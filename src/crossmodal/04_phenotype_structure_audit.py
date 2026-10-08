#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def x_number(c: str):
    m = re.fullmatch(r"X_(\d+)", c)
    return int(m.group(1)) if m else None


def phenotype_columns(path: Path):
    cols = pq.read_schema(path).names

    xs = [
        (x_number(c), c)
        for c in cols
        if x_number(c) is not None
    ]

    xs = sorted(xs)

    nums = [n for n, _ in xs]
    names = [c for _, c in xs]

    if nums:
        expected = list(
            range(min(nums), max(nums) + 1)
        )

        if nums != expected:
            raise RuntimeError(
                f"{path}: X_* columns are not contiguous."
            )

    return names


def read_dataset(
    path: Path,
    dataset: str,
    paired_ids: set[int],
):
    xcols = phenotype_columns(path)

    usecols = [
        "Metadata_NCBI_Gene_ID",
        "Metadata_Symbol",
    ] + xcols

    df = pd.read_parquet(
        path,
        columns=usecols,
    )

    df["entrez_id"] = pd.to_numeric(
        df["Metadata_NCBI_Gene_ID"],
        errors="coerce",
    ).astype("Int64")

    if df["entrez_id"].duplicated().any():
        raise RuntimeError(
            f"{dataset}: duplicate Entrez IDs."
        )

    x = df[
        df["entrez_id"].isin(paired_ids)
    ].copy()

    X = x[xcols].to_numpy(
        dtype=np.float64,
    )

    finite = np.isfinite(X)

    row_all_finite = finite.all(axis=1)

    row_norm = np.linalg.norm(
        np.nan_to_num(
            X,
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        ),
        axis=1,
    )

    row_all_zero = np.all(
        np.nan_to_num(X, nan=0.0) == 0,
        axis=1,
    )

    col_nonfinite = (
        ~finite
    ).sum(axis=0)

    col_variance = np.nanvar(
        np.where(finite, X, np.nan),
        axis=0,
    )

    summary = {
        "dataset": dataset,
        "path": str(path.resolve()),

        "n_model_rows": int(len(df)),
        "n_target_dimensions": int(len(xcols)),

        "first_target_column":
            xcols[0] if xcols else None,

        "last_target_column":
            xcols[-1] if xcols else None,

        "n_paired_rows": int(len(x)),
        "n_paired_unique_entrez":
            int(x["entrez_id"].nunique()),

        "n_rows_all_finite":
            int(row_all_finite.sum()),

        "n_rows_with_nonfinite":
            int((~row_all_finite).sum()),

        "n_all_zero_rows":
            int(row_all_zero.sum()),

        "n_columns_with_nonfinite":
            int((col_nonfinite > 0).sum()),

        "n_zero_variance_columns":
            int(
                np.sum(
                    np.isfinite(col_variance)
                    & (col_variance == 0)
                )
            ),

        "norm_summary": {
            "min": float(np.nanmin(row_norm)),
            "q01": float(
                np.nanquantile(row_norm, 0.01)
            ),
            "q25": float(
                np.nanquantile(row_norm, 0.25)
            ),
            "median": float(
                np.nanmedian(row_norm)
            ),
            "q75": float(
                np.nanquantile(row_norm, 0.75)
            ),
            "q99": float(
                np.nanquantile(row_norm, 0.99)
            ),
            "max": float(np.nanmax(row_norm)),
        },
    }

    per_gene = pd.DataFrame({
        "entrez_id": x["entrez_id"],
        "gene_symbol":
            x["Metadata_Symbol"],
        "all_finite":
            row_all_finite,
        "all_zero":
            row_all_zero,
        "l2_norm":
            row_norm,
    })

    return summary, per_gene, xcols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    root = Path(args.root).resolve()

    paired_path = (
        root
        / "results"
        / "crossmodal_v5"
        / "pairing_key_group_audit_v1"
        / "paired_genes_entrez_canonical.csv"
    )

    paired = pd.read_csv(paired_path)

    paired["entrez_id"] = pd.to_numeric(
        paired["entrez_id"],
        errors="raise",
    ).astype(int)

    paired_ids = set(
        paired["entrez_id"]
    )

    if len(paired_ids) != 5323:
        raise RuntimeError(
            f"Expected 5323 paired genes, "
            f"found {len(paired_ids)}."
        )

    outdir = (
        root
        / "results"
        / "crossmodal_v5"
        / "phenotype_structure_audit_v1"
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("[1/3] Read ORF phenotype structure")

    orf_summary, orf_gene, orf_cols = (
        read_dataset(
            root
            / "processed"
            / "orf_model_table.parquet",
            "orf",
            paired_ids,
        )
    )

    print("[2/3] Read CRISPR phenotype structure")

    crispr_summary, crispr_gene, crispr_cols = (
        read_dataset(
            root
            / "processed"
            / "crispr_model_table.parquet",
            "crispr",
            paired_ids,
        )
    )

    print("[3/3] Validate paired phenotype availability")

    orf_ids = set(
        orf_gene["entrez_id"].astype(int)
    )

    crispr_ids = set(
        crispr_gene["entrez_id"].astype(int)
    )

    missing_orf = sorted(
        paired_ids - orf_ids
    )

    missing_crispr = sorted(
        paired_ids - crispr_ids
    )

    overlap = orf_gene.merge(
        crispr_gene,
        on="entrez_id",
        how="inner",
        suffixes=("_orf", "_crispr"),
        validate="one_to_one",
    )

    overlap.to_csv(
        outdir
        / "paired_phenotype_qc.csv",
        index=False,
    )

    status = (
        "PASS"
        if (
            len(missing_orf) == 0
            and len(missing_crispr) == 0
            and orf_summary[
                "n_rows_with_nonfinite"
            ] == 0
            and crispr_summary[
                "n_rows_with_nonfinite"
            ] == 0
        )
        else "REVIEW"
    )

    summary = {
        "status":
            f"PHENOTYPE_STRUCTURE_{status}",

        "split_protocol_already_frozen":
            True,

        "split_changes_allowed":
            False,

        "performance_metric_computed":
            False,

        "prediction_model_fit":
            False,

        "paired_universe_n":
            5323,

        "paired_phenotype_overlap_n":
            int(len(overlap)),

        "missing_orf_n":
            len(missing_orf),

        "missing_crispr_n":
            len(missing_crispr),

        "missing_orf_entrez":
            missing_orf,

        "missing_crispr_entrez":
            missing_crispr,

        "orf":
            orf_summary,

        "crispr":
            crispr_summary,

        "target_space_note": (
            "ORF and CRISPR target spaces are "
            "kept modality-specific. Different "
            "dimensions do not prevent supervised "
            "cross-space prediction."
        ),

        "important": (
            "No train/test performance, correlation, "
            "R2, cosine prediction metric, or model "
            "fitting was performed."
        ),
    }

    with open(
        outdir / "summary.json",
        "w",
    ) as f:
        json.dump(
            summary,
            f,
            indent=2,
            ensure_ascii=False,
        )

    with open(
        outdir / "orf_target_columns.txt",
        "w",
    ) as f:
        f.write("\n".join(orf_cols))

    with open(
        outdir / "crispr_target_columns.txt",
        "w",
    ) as f:
        f.write("\n".join(crispr_cols))

    print()
    print(
        json.dumps(
            summary,
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("[DONE]", outdir)


if __name__ == "__main__":
    main()
