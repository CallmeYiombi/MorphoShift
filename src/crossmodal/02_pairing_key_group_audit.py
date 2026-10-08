#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


BENCHMARKS = {
    "random_gene": None,
    "homology30": "cluster_id30",
    "homology50": "cluster_id50",
    "panther_family": "panther_family",
}

SEEDS = [42, 123, 2024]


def norm_symbol(x):
    x = x.astype("string").str.strip().str.upper()
    return x.mask(x.isin(["", "NAN", "NONE", "<NA>"]))


def norm_uniprot(x):
    x = x.astype("string").str.strip().str.upper()
    x = x.mask(x.isin(["", "NAN", "NONE", "<NA>"]))
    return x


def read_model_identity(path: Path, dataset: str):
    cols = [
        "Metadata_NCBI_Gene_ID",
        "Metadata_Symbol",
        "uniprot_accession",
    ]

    df = pd.read_parquet(path, columns=cols)

    out = pd.DataFrame({
        "entrez_id": pd.to_numeric(
            df["Metadata_NCBI_Gene_ID"],
            errors="coerce",
        ).astype("Int64"),
        "gene_symbol": norm_symbol(df["Metadata_Symbol"]),
        "uniprot_accession": norm_uniprot(df["uniprot_accession"]),
    })

    out["uniprot_canonical"] = (
        out["uniprot_accession"]
        .str.split("-", n=1)
        .str[0]
    )

    if out["entrez_id"].isna().any():
        raise RuntimeError(
            f"{dataset}: missing Entrez IDs found."
        )

    if out["entrez_id"].duplicated().any():
        dup = out.loc[
            out["entrez_id"].duplicated(False)
        ].sort_values("entrez_id")
        raise RuntimeError(
            f"{dataset}: duplicate Entrez IDs found:\n"
            f"{dup.head(20).to_string(index=False)}"
        )

    return out


def read_split(path: Path, group_col: str | None):
    df = pd.read_csv(path)

    required = [
        "Metadata_NCBI_Gene_ID",
        "Metadata_Symbol",
        "split",
    ]

    if group_col is not None:
        required.append(group_col)

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"{path}: missing columns {missing}"
        )

    out = pd.DataFrame({
        "entrez_id": pd.to_numeric(
            df["Metadata_NCBI_Gene_ID"],
            errors="coerce",
        ).astype("Int64"),
        "gene_symbol": norm_symbol(df["Metadata_Symbol"]),
        "split": df["split"].astype("string"),
    })

    if group_col is not None:
        out["group"] = df[group_col].astype("string")
        out.loc[
            out["group"].isin(["", "nan", "None", "<NA>"]),
            "group"
        ] = pd.NA

    return out


def add_group_signatures(df: pd.DataFrame):
    """
    Compare partitions without assuming ORF/CRISPR group labels
    themselves use identical numbering.

    Each gene gets the sorted set of paired Entrez IDs belonging to
    its group.
    """
    x = df.dropna(subset=["group"]).copy()

    if len(x) == 0:
        return pd.DataFrame(
            columns=["entrez_id", "group", "group_signature"]
        )

    members = (
        x.groupby("group")["entrez_id"]
        .apply(
            lambda s: "|".join(
                map(
                    str,
                    sorted(
                        int(v)
                        for v in s.dropna().unique()
                    ),
                )
            )
        )
        .rename("group_signature")
    )

    return x.merge(
        members,
        left_on="group",
        right_index=True,
        how="left",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    root = Path(args.root).resolve()

    outdir = (
        root
        / "results"
        / "crossmodal_v5"
        / "pairing_key_group_audit_v1"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Load model identities using Entrez as canonical key")

    orf = read_model_identity(
        root / "processed/orf_model_table.parquet",
        "orf",
    )

    crispr = read_model_identity(
        root / "processed/crispr_model_table.parquet",
        "crispr",
    )

    paired = orf.merge(
        crispr,
        on="entrez_id",
        how="inner",
        suffixes=("_orf", "_crispr"),
        validate="one_to_one",
    )

    paired["symbol_match"] = (
        paired["gene_symbol_orf"]
        == paired["gene_symbol_crispr"]
    )

    paired["uniprot_exact_match"] = (
        paired["uniprot_accession_orf"].notna()
        & paired["uniprot_accession_crispr"].notna()
        & (
            paired["uniprot_accession_orf"]
            == paired["uniprot_accession_crispr"]
        )
    )

    paired["uniprot_canonical_match"] = (
        paired["uniprot_canonical_orf"].notna()
        & paired["uniprot_canonical_crispr"].notna()
        & (
            paired["uniprot_canonical_orf"]
            == paired["uniprot_canonical_crispr"]
        )
    )

    paired.to_csv(
        outdir / "paired_genes_entrez_canonical.csv",
        index=False,
    )

    paired.loc[~paired["symbol_match"]].to_csv(
        outdir / "paired_symbol_aliases.csv",
        index=False,
    )

    paired.loc[
        ~paired["uniprot_canonical_match"]
    ].to_csv(
        outdir / "paired_uniprot_disagreements.csv",
        index=False,
    )

    print(
        f"  paired Entrez genes={len(paired)}"
    )
    print(
        f"  symbol aliases="
        f"{int((~paired['symbol_match']).sum())}"
    )
    print(
        f"  canonical UniProt disagreements="
        f"{int((~paired['uniprot_canonical_match']).sum())}"
    )

    paired_ids = set(
        paired["entrez_id"].astype(int)
    )

    print("[2/5] Verify group definitions are seed-stable")

    seed_stability_records = []

    for dataset in ["orf", "crispr"]:
        for benchmark, group_col in BENCHMARKS.items():

            if group_col is None:
                continue

            tables = {}

            for seed in SEEDS:
                path = (
                    root
                    / "processed"
                    / "splits"
                    / dataset
                    / benchmark
                    / f"seed_{seed}.csv"
                )

                x = read_split(path, group_col)
                x = x[
                    x["entrez_id"].isin(paired_ids)
                ].copy()

                tables[seed] = x[
                    ["entrez_id", "group"]
                ].rename(
                    columns={"group": f"group_{seed}"}
                )

            z = tables[SEEDS[0]]

            for seed in SEEDS[1:]:
                z = z.merge(
                    tables[seed],
                    on="entrez_id",
                    how="outer",
                    validate="one_to_one",
                )

            base = f"group_{SEEDS[0]}"

            comparable = z[base].notna()

            for seed in SEEDS[1:]:
                comparable &= z[f"group_{seed}"].notna()

            zz = z.loc[comparable].copy()

            if len(zz):
                stable = np.ones(len(zz), dtype=bool)

                for seed in SEEDS[1:]:
                    stable &= (
                        zz[base].astype(str).to_numpy()
                        == zz[f"group_{seed}"]
                        .astype(str)
                        .to_numpy()
                    )

                n_stable = int(stable.sum())
                frac = float(stable.mean())
            else:
                n_stable = 0
                frac = None

            seed_stability_records.append({
                "dataset": dataset,
                "benchmark": benchmark,
                "n_paired_total": len(paired),
                "n_comparable": len(zz),
                "n_seed_stable": n_stable,
                "seed_stability_fraction": frac,
            })

    seed_stability = pd.DataFrame(
        seed_stability_records
    )

    seed_stability.to_csv(
        outdir / "group_seed_stability.csv",
        index=False,
    )

    print("[3/5] Compare ORF vs CRISPR biological partitions")

    partition_records = []

    # Group definitions should be seed-independent.
    # Use seed 42 after explicitly auditing stability above.
    for benchmark, group_col in BENCHMARKS.items():

        if group_col is None:
            continue

        opath = (
            root
            / "processed"
            / "splits"
            / "orf"
            / benchmark
            / "seed_42.csv"
        )

        cpath = (
            root
            / "processed"
            / "splits"
            / "crispr"
            / benchmark
            / "seed_42.csv"
        )

        o = read_split(opath, group_col)
        c = read_split(cpath, group_col)

        o = o[
            o["entrez_id"].isin(paired_ids)
        ].copy()

        c = c[
            c["entrez_id"].isin(paired_ids)
        ].copy()

        common_ids = (
            set(o.loc[o["group"].notna(), "entrez_id"])
            & set(c.loc[c["group"].notna(), "entrez_id"])
        )

        o2 = o[
            o["entrez_id"].isin(common_ids)
        ][["entrez_id", "group"]].copy()

        c2 = c[
            c["entrez_id"].isin(common_ids)
        ][["entrez_id", "group"]].copy()

        osig = add_group_signatures(o2)[
            ["entrez_id", "group", "group_signature"]
        ].rename(
            columns={
                "group": "group_orf",
                "group_signature": "signature_orf",
            }
        )

        csig = add_group_signatures(c2)[
            ["entrez_id", "group", "group_signature"]
        ].rename(
            columns={
                "group": "group_crispr",
                "group_signature": "signature_crispr",
            }
        )

        z = osig.merge(
            csig,
            on="entrez_id",
            how="inner",
            validate="one_to_one",
        )

        z["partition_match"] = (
            z["signature_orf"]
            == z["signature_crispr"]
        )

        z.to_csv(
            outdir
            / f"{benchmark}_crossmodal_partition_gene_audit.csv",
            index=False,
        )

        partition_records.append({
            "benchmark": benchmark,
            "n_paired_total": len(paired),
            "n_orf_with_group": int(
                o["group"].notna().sum()
            ),
            "n_crispr_with_group": int(
                c["group"].notna().sum()
            ),
            "n_common_group_eligible": len(z),
            "n_partition_match": int(
                z["partition_match"].sum()
            ),
            "partition_match_fraction":
                float(z["partition_match"].mean())
                if len(z)
                else None,
        })

    partition_df = pd.DataFrame(
        partition_records
    )

    partition_df.to_csv(
        outdir / "crossmodal_group_partition_concordance.csv",
        index=False,
    )

    print("[4/5] Audit current ORF/CRISPR split concordance")

    split_records = []

    for benchmark, group_col in BENCHMARKS.items():

        for seed in SEEDS:

            opath = (
                root
                / "processed"
                / "splits"
                / "orf"
                / benchmark
                / f"seed_{seed}.csv"
            )

            cpath = (
                root
                / "processed"
                / "splits"
                / "crispr"
                / benchmark
                / f"seed_{seed}.csv"
            )

            o = read_split(opath, group_col)
            c = read_split(cpath, group_col)

            o = o[
                o["entrez_id"].isin(paired_ids)
            ][["entrez_id", "split"]].rename(
                columns={"split": "split_orf"}
            )

            c = c[
                c["entrez_id"].isin(paired_ids)
            ][["entrez_id", "split"]].rename(
                columns={"split": "split_crispr"}
            )

            z = o.merge(
                c,
                on="entrez_id",
                how="inner",
                validate="one_to_one",
            )

            z["same_split"] = (
                z["split_orf"]
                == z["split_crispr"]
            )

            split_records.append({
                "benchmark": benchmark,
                "seed": seed,
                "n_paired_total": len(paired),
                "n_comparable": len(z),
                "same_split_n":
                    int(z["same_split"].sum()),
                "same_split_fraction":
                    float(z["same_split"].mean())
                    if len(z)
                    else None,
            })

    split_df = pd.DataFrame(split_records)

    split_df.to_csv(
        outdir / "existing_split_crossmodal_concordance.csv",
        index=False,
    )

    print("[5/5] Write summary")

    summary = {
        "status": "PAIRING_KEY_GROUP_AUDIT_COMPLETE",
        "performance_unblinded": False,
        "phenotype_columns_read": False,

        "canonical_pairing_key":
            "Metadata_NCBI_Gene_ID",

        "n_paired_entrez": int(len(paired)),

        "n_symbol_match": int(
            paired["symbol_match"].sum()
        ),

        "n_symbol_alias": int(
            (~paired["symbol_match"]).sum()
        ),

        "n_uniprot_exact_match": int(
            paired["uniprot_exact_match"].sum()
        ),

        "n_uniprot_canonical_match": int(
            paired["uniprot_canonical_match"].sum()
        ),

        "fraction_of_orf": float(
            len(paired) / len(orf)
        ),

        "fraction_of_crispr": float(
            len(paired) / len(crispr)
        ),

        "group_seed_stability":
            seed_stability.to_dict(orient="records"),

        "crossmodal_partition_concordance":
            partition_df.to_dict(orient="records"),

        "existing_split_crossmodal_concordance":
            split_df.to_dict(orient="records"),

        "recommendation": (
            "Do not reuse modality-specific split assignments "
            "for MorphoShift v5. After group concordance is "
            "validated, generate a new paired-gene split with "
            "identical train/validation/test genes for both "
            "ORF->CRISPR and CRISPR->ORF directions."
        ),
    }

    with open(outdir / "summary.json", "w") as f:
        json.dump(
            summary,
            f,
            indent=2,
            ensure_ascii=False,
        )

    print()
    print("=== CANONICAL PAIRING ===")
    print(f"paired genes             : {len(paired)}")
    print(
        "symbol aliases           :",
        int((~paired["symbol_match"]).sum()),
    )
    print(
        "canonical UniProt mismatch:",
        int(
            (~paired["uniprot_canonical_match"]).sum()
        ),
    )

    print()
    print("=== GROUP SEED STABILITY ===")
    print(seed_stability.to_string(index=False))

    print()
    print("=== CROSS-MODAL PARTITION CONCORDANCE ===")
    print(partition_df.to_string(index=False))

    print()
    print("=== EXISTING SPLIT CONCORDANCE ===")
    print(split_df.to_string(index=False))

    print()
    print("[DONE]", outdir)
    print(
        "[IMPORTANT] No morphology X_* column was read."
    )


if __name__ == "__main__":
    main()
