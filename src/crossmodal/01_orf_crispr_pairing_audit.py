#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def choose_column(columns, candidates):
    lower = {c.lower(): c for c in columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def norm_symbol(s):
    s = s.astype("string").str.strip().str.upper()
    return s.mask(s.isin(["", "NAN", "NONE", "<NA>"]))


def norm_uniprot(s):
    s = s.astype("string").str.strip().str.upper()
    return s.mask(s.isin(["", "NAN", "NONE", "<NA>"]))


def read_identity(path: Path, dataset: str):
    cols = pq.read_schema(path).names

    gene_col = choose_column(
        cols,
        [
            "Metadata_Symbol",
            "gene_symbol",
            "symbol",
            "gene",
            "Gene",
        ],
    )

    entrez_col = choose_column(
        cols,
        [
            "Metadata_NCBI_Gene_ID",
            "ncbi_gene_id",
            "entrez_id",
            "gene_id",
        ],
    )

    uniprot_col = choose_column(
        cols,
        [
            "uniprot_accession",
            "UniProt",
            "uniprot",
            "protein_accession",
        ],
    )

    prot_match_col = choose_column(
        cols,
        [
            "Metadata_Prot_Match",
            "Prot_Match",
            "orf_prot_match",
        ],
    )

    if gene_col is None:
        raise RuntimeError(
            f"{dataset}: gene-symbol column not found in {path}"
        )

    usecols = [gene_col]

    for c in [entrez_col, uniprot_col, prot_match_col]:
        if c is not None and c not in usecols:
            usecols.append(c)

    # IMPORTANT:
    # Only identity columns are loaded.
    # X_1 ... morphology columns are never read here.
    df = pd.read_parquet(path, columns=usecols)

    out = pd.DataFrame()
    out["gene_symbol"] = norm_symbol(df[gene_col])

    if entrez_col is not None:
        out["entrez_id"] = pd.to_numeric(
            df[entrez_col], errors="coerce"
        ).astype("Int64")
    else:
        out["entrez_id"] = pd.Series(
            pd.NA, index=df.index, dtype="Int64"
        )

    if uniprot_col is not None:
        out["uniprot_accession"] = norm_uniprot(df[uniprot_col])
        out["uniprot_canonical"] = (
            out["uniprot_accession"]
            .str.split("-", n=1)
            .str[0]
        )
    else:
        out["uniprot_accession"] = pd.Series(
            pd.NA, index=df.index, dtype="string"
        )
        out["uniprot_canonical"] = pd.Series(
            pd.NA, index=df.index, dtype="string"
        )

    if prot_match_col is not None:
        out["prot_match"] = df[prot_match_col]
    else:
        out["prot_match"] = pd.NA

    meta = {
        "dataset": dataset,
        "path": str(path.resolve()),
        "n_rows": int(len(out)),
        "n_gene_nonnull": int(out["gene_symbol"].notna().sum()),
        "n_gene_unique": int(out["gene_symbol"].nunique(dropna=True)),
        "n_entrez_nonnull": int(out["entrez_id"].notna().sum()),
        "n_entrez_unique": int(out["entrez_id"].nunique(dropna=True)),
        "n_uniprot_nonnull": int(
            out["uniprot_accession"].notna().sum()
        ),
        "n_uniprot_unique": int(
            out["uniprot_accession"].nunique(dropna=True)
        ),
        "gene_col": gene_col,
        "entrez_col": entrez_col,
        "uniprot_col": uniprot_col,
        "prot_match_col": prot_match_col,
        "phenotype_columns_read": False,
    }

    return out, meta


def set_overlap(a, b, col):
    aa = set(a[col].dropna().tolist())
    bb = set(b[col].dropna().tolist())
    return len(aa & bb)


def inspect_split_files(root: Path):
    records = []
    details = {}

    base = root / "processed" / "splits"

    for dataset in ["orf", "crispr"]:
        ds_root = base / dataset

        if not ds_root.exists():
            continue

        for path in sorted(ds_root.rglob("seed_*.csv")):
            rel = path.relative_to(root)
            df = pd.read_csv(path)

            rel_ds = path.relative_to(ds_root)
            benchmark = (
                rel_ds.parts[0]
                if len(rel_ds.parts) > 1
                else "UNKNOWN"
            )

            m = re.search(r"seed_(\d+)", path.stem)
            seed = int(m.group(1)) if m else None

            candidate_split_cols = [
                c for c in df.columns
                if c.lower() in {
                    "split",
                    "set",
                    "partition",
                    "subset",
                    "fold",
                }
            ]

            records.append(
                {
                    "dataset": dataset,
                    "benchmark": benchmark,
                    "seed": seed,
                    "path": str(rel),
                    "n_rows": len(df),
                    "n_columns": len(df.columns),
                    "columns": "|".join(map(str, df.columns)),
                    "candidate_split_columns":
                        "|".join(candidate_split_cols),
                }
            )

            details[str(rel)] = {
                "columns": list(map(str, df.columns)),
                "candidate_split_value_counts": {
                    c: {
                        str(k): int(v)
                        for k, v in
                        df[c]
                        .astype("string")
                        .value_counts(dropna=False)
                        .to_dict()
                        .items()
                    }
                    for c in candidate_split_cols
                },
            }

    return pd.DataFrame(records), details


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    root = Path(args.root).resolve()

    orf_path = root / "processed" / "orf_model_table.parquet"
    crispr_path = root / "processed" / "crispr_model_table.parquet"

    outdir = (
        root
        / "results"
        / "crossmodal_v5"
        / "blind_pairing_audit_v1"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Load identity columns only")

    orf, orf_meta = read_identity(orf_path, "orf")
    crispr, crispr_meta = read_identity(crispr_path, "crispr")

    print(
        f"  ORF:    rows={len(orf)} "
        f"genes={orf.gene_symbol.nunique()}"
    )
    print(
        f"  CRISPR: rows={len(crispr)} "
        f"genes={crispr.gene_symbol.nunique()}"
    )

    if orf["gene_symbol"].duplicated().any():
        raise RuntimeError(
            "ORF contains duplicate normalized gene symbols."
        )

    if crispr["gene_symbol"].duplicated().any():
        raise RuntimeError(
            "CRISPR contains duplicate normalized gene symbols."
        )

    print("[2/5] Build exact gene-symbol paired universe")

    paired = orf.merge(
        crispr,
        on="gene_symbol",
        how="inner",
        suffixes=("_orf", "_crispr"),
        validate="one_to_one",
    )

    paired["entrez_both"] = (
        paired["entrez_id_orf"].notna()
        & paired["entrez_id_crispr"].notna()
    )

    paired["entrez_match"] = (
        paired["entrez_both"]
        & (
            paired["entrez_id_orf"]
            == paired["entrez_id_crispr"]
        )
    )

    paired["uniprot_both"] = (
        paired["uniprot_accession_orf"].notna()
        & paired["uniprot_accession_crispr"].notna()
    )

    paired["uniprot_exact_match"] = (
        paired["uniprot_both"]
        & (
            paired["uniprot_accession_orf"]
            == paired["uniprot_accession_crispr"]
        )
    )

    paired["uniprot_canonical_both"] = (
        paired["uniprot_canonical_orf"].notna()
        & paired["uniprot_canonical_crispr"].notna()
    )

    paired["uniprot_canonical_match"] = (
        paired["uniprot_canonical_both"]
        & (
            paired["uniprot_canonical_orf"]
            == paired["uniprot_canonical_crispr"]
        )
    )

    print("[3/5] Audit cross-key disagreements")

    symbol_conflicts = paired[
        (
            paired["entrez_both"]
            & ~paired["entrez_match"]
        )
        |
        (
            paired["uniprot_canonical_both"]
            & ~paired["uniprot_canonical_match"]
        )
    ].copy()

    # Same Entrez ID but different current symbols
    entrez_cross = (
        orf[
            ["gene_symbol", "entrez_id"]
        ]
        .dropna()
        .merge(
            crispr[
                ["gene_symbol", "entrez_id"]
            ].dropna(),
            on="entrez_id",
            how="inner",
            suffixes=("_orf", "_crispr"),
        )
    )

    entrez_symbol_mismatch = entrez_cross[
        entrez_cross["gene_symbol_orf"]
        != entrez_cross["gene_symbol_crispr"]
    ].copy()

    # Same canonical UniProt but different symbols
    up_cross = (
        orf[
            ["gene_symbol", "uniprot_canonical"]
        ]
        .dropna()
        .merge(
            crispr[
                ["gene_symbol", "uniprot_canonical"]
            ].dropna(),
            on="uniprot_canonical",
            how="inner",
            suffixes=("_orf", "_crispr"),
        )
    )

    uniprot_symbol_mismatch = up_cross[
        up_cross["gene_symbol_orf"]
        != up_cross["gene_symbol_crispr"]
    ].copy()

    print("[4/5] Inventory canonical frozen split files")

    split_inventory, split_details = inspect_split_files(root)

    print("[5/5] Write audit outputs")

    n_pair = len(paired)

    summary = {
        "status": "BLIND_PAIRING_AUDIT_COMPLETE",
        "performance_unblinded": False,
        "phenotype_columns_read": False,
        "primary_pairing_key": "normalized exact gene symbol",
        "orf": orf_meta,
        "crispr": crispr_meta,
        "pairing": {
            "n_symbol_pairs": int(n_pair),
            "fraction_of_orf":
                float(n_pair / len(orf)) if len(orf) else None,
            "fraction_of_crispr":
                float(n_pair / len(crispr)) if len(crispr) else None,

            "n_entrez_overlap_independent":
                int(set_overlap(orf, crispr, "entrez_id")),

            "n_uniprot_exact_overlap_independent":
                int(
                    set_overlap(
                        orf,
                        crispr,
                        "uniprot_accession",
                    )
                ),

            "n_uniprot_canonical_overlap_independent":
                int(
                    set_overlap(
                        orf,
                        crispr,
                        "uniprot_canonical",
                    )
                ),

            "n_symbol_pairs_both_entrez":
                int(paired["entrez_both"].sum()),

            "n_symbol_pairs_entrez_match":
                int(paired["entrez_match"].sum()),

            "n_symbol_pairs_entrez_mismatch":
                int(
                    (
                        paired["entrez_both"]
                        & ~paired["entrez_match"]
                    ).sum()
                ),

            "n_symbol_pairs_both_uniprot":
                int(paired["uniprot_both"].sum()),

            "n_symbol_pairs_uniprot_exact_match":
                int(paired["uniprot_exact_match"].sum()),

            "n_symbol_pairs_uniprot_canonical_match":
                int(
                    paired[
                        "uniprot_canonical_match"
                    ].sum()
                ),

            "n_symbol_pair_identity_conflicts":
                int(len(symbol_conflicts)),

            "n_same_entrez_different_symbol":
                int(len(entrez_symbol_mismatch)),

            "n_same_uniprot_different_symbol":
                int(len(uniprot_symbol_mismatch)),
        },
        "split_policy_note": (
            "Existing ORF and CRISPR split artifacts are only "
            "inventoried here. No split semantics are inferred and "
            "no cross-modal train/test split is created in this audit."
        ),
        "next_gate": (
            "Define and freeze the paired-gene MorphoShift v5 "
            "split protocol before reading morphology target columns."
        ),
    }

    paired.to_csv(
        outdir / "paired_genes_identity.csv",
        index=False,
    )

    symbol_conflicts.to_csv(
        outdir / "symbol_pair_identity_conflicts.csv",
        index=False,
    )

    entrez_symbol_mismatch.to_csv(
        outdir / "same_entrez_different_symbol.csv",
        index=False,
    )

    uniprot_symbol_mismatch.to_csv(
        outdir / "same_uniprot_different_symbol.csv",
        index=False,
    )

    split_inventory.to_csv(
        outdir / "canonical_split_inventory.csv",
        index=False,
    )

    with open(
        outdir / "canonical_split_inventory_details.json",
        "w",
    ) as f:
        json.dump(
            split_details,
            f,
            indent=2,
            ensure_ascii=False,
        )

    with open(outdir / "summary.json", "w") as f:
        json.dump(
            summary,
            f,
            indent=2,
            ensure_ascii=False,
        )

    print()
    print("=== BLIND CROSS-MODAL PAIRING SUMMARY ===")
    print(
        json.dumps(
            summary["pairing"],
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("[DONE]", outdir)
    print(
        "[IMPORTANT] No X_* morphology target column "
        "was read."
    )


if __name__ == "__main__":
    main()
