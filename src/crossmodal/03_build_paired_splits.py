#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


SEEDS = [42, 123, 2024]

RATIOS = {
    "train": 0.80,
    "val": 0.10,
    "test": 0.10,
}

SPLITS = ["train", "val", "test"]

BENCHMARKS = {
    "random_gene": None,
    "homology30": "cluster_id30",
    "homology50": "cluster_id50",
    "panther_family": "panther_family",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def sha256_split(df: pd.DataFrame) -> str:
    x = (
        df[["entrez_id", "split"]]
        .sort_values("entrez_id")
        .astype(str)
    )

    payload = "\n".join(
        x["entrez_id"] + "\t" + x["split"]
    ).encode()

    return hashlib.sha256(payload).hexdigest()


def integer_targets(n: int) -> dict[str, int]:
    """
    Largest-remainder allocation for 80/10/10.
    """
    exact = {
        k: n * RATIOS[k]
        for k in SPLITS
    }

    base = {
        k: int(np.floor(exact[k]))
        for k in SPLITS
    }

    remainder = n - sum(base.values())

    ranked = sorted(
        SPLITS,
        key=lambda k: exact[k] - base[k],
        reverse=True,
    )

    for k in ranked[:remainder]:
        base[k] += 1

    assert sum(base.values()) == n
    return base


def make_random_gene_split(
    genes: pd.DataFrame,
    seed: int,
) -> pd.DataFrame:

    genes = genes.copy().reset_index(drop=True)

    if genes["entrez_id"].duplicated().any():
        raise RuntimeError(
            "Duplicate Entrez IDs in paired universe."
        )

    rng = np.random.default_rng(seed)

    perm = rng.permutation(len(genes))

    targets = integer_targets(len(genes))

    labels = np.empty(len(genes), dtype=object)

    start = 0

    for split in SPLITS:
        n = targets[split]
        idx = perm[start:start + n]
        labels[idx] = split
        start += n

    assert start == len(genes)

    genes["split"] = labels

    return genes


def grouped_assignment(
    genes: pd.DataFrame,
    group_col: str,
    seed: int,
    n_restarts: int = 256,
) -> pd.DataFrame:
    """
    Assign whole biological groups to train/val/test.

    Optimization target:
      approximate 80/10/10 by number of genes,
      while never splitting a group.

    Multiple deterministic randomized restarts are used and
    the allocation with the lowest normalized count error is kept.
    """

    if genes[group_col].isna().any():
        raise RuntimeError(
            f"Missing {group_col} values present."
        )

    group_sizes = (
        genes.groupby(group_col, sort=False)
        .size()
        .rename("n")
        .reset_index()
    )

    target = {
        k: RATIOS[k] * len(genes)
        for k in SPLITS
    }

    best_assignment = None
    best_score = None

    base_rng = np.random.default_rng(seed)

    restart_seeds = base_rng.integers(
        0,
        2**32 - 1,
        size=n_restarts,
        dtype=np.uint32,
    )

    for restart_seed in restart_seeds:

        rng = np.random.default_rng(
            int(restart_seed)
        )

        work = group_sizes.copy()

        # Random jitter gives seed-dependent solutions,
        # while sorting large groups first improves balance.
        work["_jitter"] = rng.random(len(work))

        work = work.sort_values(
            ["n", "_jitter"],
            ascending=[False, True],
        ).reset_index(drop=True)

        counts = {
            k: 0
            for k in SPLITS
        }

        assignment = {}

        for _, row in work.iterrows():

            group = row[group_col]
            size = int(row["n"])

            candidate_scores = []

            split_order = list(
                rng.permutation(SPLITS)
            )

            for split in split_order:

                proposed = counts.copy()
                proposed[split] += size

                score = sum(
                    (
                        (
                            proposed[k]
                            - target[k]
                        )
                        / max(target[k], 1.0)
                    ) ** 2
                    for k in SPLITS
                )

                # Mild extra penalty for overshooting
                # validation/test targets.
                for k in ["val", "test"]:
                    if proposed[k] > target[k]:
                        overshoot = (
                            proposed[k]
                            - target[k]
                        ) / max(target[k], 1.0)

                        score += (
                            0.25
                            * overshoot**2
                        )

                candidate_scores.append(
                    (score, split)
                )

            score, chosen = min(
                candidate_scores,
                key=lambda x: x[0],
            )

            assignment[group] = chosen
            counts[chosen] += size

        final_score = sum(
            (
                (
                    counts[k]
                    - target[k]
                )
                / max(target[k], 1.0)
            ) ** 2
            for k in SPLITS
        )

        if (
            best_score is None
            or final_score < best_score
        ):
            best_score = final_score
            best_assignment = assignment

    out = genes.copy()

    out["split"] = (
        out[group_col]
        .map(best_assignment)
    )

    if out["split"].isna().any():
        raise RuntimeError(
            "Some groups were not assigned."
        )

    return out


def read_group_map(
    root: Path,
    benchmark: str,
    group_col: str,
) -> pd.DataFrame:
    """
    Biological group definitions were already audited as:
      - seed stable
      - ORF/CRISPR partition concordant

    Use the ORF seed-42 artifact only as the canonical
    source for the group label itself.
    """

    path = (
        root
        / "processed"
        / "splits"
        / "orf"
        / benchmark
        / "seed_42.csv"
    )

    df = pd.read_csv(path)

    required = {
        "Metadata_NCBI_Gene_ID",
        group_col,
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"{path}: missing {sorted(missing)}"
        )

    out = pd.DataFrame({
        "entrez_id": pd.to_numeric(
            df["Metadata_NCBI_Gene_ID"],
            errors="coerce",
        ).astype("Int64"),

        group_col:
            df[group_col].astype("string"),
    })

    out = out.dropna(
        subset=["entrez_id", group_col]
    )

    if out["entrez_id"].duplicated().any():
        raise RuntimeError(
            f"{benchmark}: duplicate Entrez IDs."
        )

    return out


def audit_split(
    df: pd.DataFrame,
    benchmark: str,
    seed: int,
    group_col: str | None,
) -> dict:

    counts = (
        df["split"]
        .value_counts()
        .reindex(SPLITS, fill_value=0)
    )

    record = {
        "benchmark": benchmark,
        "seed": seed,
        "n_genes": int(len(df)),
        "n_train": int(counts["train"]),
        "n_val": int(counts["val"]),
        "n_test": int(counts["test"]),
        "frac_train":
            float(counts["train"] / len(df)),
        "frac_val":
            float(counts["val"] / len(df)),
        "frac_test":
            float(counts["test"] / len(df)),
        "split_sha256":
            sha256_split(df),
    }

    if group_col is None:
        record.update({
            "n_groups": int(len(df)),
            "group_leakage_count": 0,
        })

    else:
        group_split_counts = (
            df.groupby(group_col)["split"]
            .nunique()
        )

        leakage = (
            group_split_counts > 1
        ).sum()

        record.update({
            "n_groups":
                int(df[group_col].nunique()),
            "group_leakage_count":
                int(leakage),
        })

    return record


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        default=".",
    )

    ap.add_argument(
        "--n-restarts",
        type=int,
        default=256,
    )

    args = ap.parse_args()

    root = Path(args.root).resolve()

    pair_path = (
        root
        / "results"
        / "crossmodal_v5"
        / "pairing_key_group_audit_v1"
        / "paired_genes_entrez_canonical.csv"
    )

    if not pair_path.exists():
        raise RuntimeError(
            f"Missing paired universe: {pair_path}"
        )

    print("[1/5] Load frozen paired Entrez universe")

    paired = pd.read_csv(pair_path)

    paired["entrez_id"] = pd.to_numeric(
        paired["entrez_id"],
        errors="raise",
    ).astype("Int64")

    if len(paired) != 5323:
        raise RuntimeError(
            f"Expected 5323 paired genes, "
            f"found {len(paired)}."
        )

    if paired["entrez_id"].duplicated().any():
        raise RuntimeError(
            "Duplicate paired Entrez IDs."
        )

    base_cols = [
        "entrez_id",
        "gene_symbol_orf",
        "gene_symbol_crispr",
        "uniprot_accession_orf",
        "uniprot_accession_crispr",
    ]

    missing = [
        c for c in base_cols
        if c not in paired.columns
    ]

    if missing:
        raise RuntimeError(
            f"Paired file missing columns: "
            f"{missing}"
        )

    paired_base = paired[
        base_cols
    ].copy()

    outroot = (
        root
        / "processed"
        / "splits"
        / "crossmodal_v5"
    )

    outroot.mkdir(
        parents=True,
        exist_ok=True,
    )

    audit_records = []
    source_manifest = {}

    print("[2/5] Load canonical biological groups")

    group_maps = {}

    for benchmark, group_col in (
        BENCHMARKS.items()
    ):

        if group_col is None:
            continue

        group_map = read_group_map(
            root,
            benchmark,
            group_col,
        )

        group_maps[benchmark] = group_map

        source_path = (
            root
            / "processed"
            / "splits"
            / "orf"
            / benchmark
            / "seed_42.csv"
        )

        source_manifest[benchmark] = {
            "path":
                str(source_path.relative_to(root)),
            "sha256":
                sha256_file(source_path),
            "group_column":
                group_col,
        }

    print("[3/5] Generate paired-gene splits")

    for benchmark, group_col in (
        BENCHMARKS.items()
    ):

        benchmark_dir = (
            outroot / benchmark
        )

        benchmark_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        if group_col is None:

            universe = paired_base.copy()

        else:

            universe = paired_base.merge(
                group_maps[benchmark],
                on="entrez_id",
                how="inner",
                validate="one_to_one",
            )

            if benchmark in {
                "homology30",
                "homology50",
            }:
                if len(universe) != 5323:
                    raise RuntimeError(
                        f"{benchmark}: expected "
                        f"5323 eligible genes, "
                        f"found {len(universe)}."
                    )

            if benchmark == "panther_family":
                if len(universe) != 5303:
                    raise RuntimeError(
                        "panther_family: expected "
                        f"5303 eligible genes, "
                        f"found {len(universe)}."
                    )

        print(
            f"  {benchmark}: "
            f"eligible genes={len(universe)}"
        )

        for seed in SEEDS:

            if group_col is None:

                split_df = (
                    make_random_gene_split(
                        universe,
                        seed,
                    )
                )

            else:

                split_df = grouped_assignment(
                    universe,
                    group_col,
                    seed,
                    n_restarts=args.n_restarts,
                )

            # Stable output ordering
            split_df = (
                split_df
                .sort_values("entrez_id")
                .reset_index(drop=True)
            )

            split_df.insert(
                0,
                "split_protocol",
                "MorphoShift_v5_paired",
            )

            split_df.insert(
                1,
                "benchmark",
                benchmark,
            )

            split_df.insert(
                2,
                "seed",
                seed,
            )

            outpath = (
                benchmark_dir
                / f"seed_{seed}.csv"
            )

            split_df.to_csv(
                outpath,
                index=False,
            )

            audit = audit_split(
                split_df,
                benchmark,
                seed,
                group_col,
            )

            audit["path"] = str(
                outpath.relative_to(root)
            )

            audit_records.append(audit)

    print("[4/5] Validate all split artifacts")

    audit_df = pd.DataFrame(
        audit_records
    )

    if (
        audit_df[
            "group_leakage_count"
        ] != 0
    ).any():
        raise RuntimeError(
            "Group leakage detected."
        )

    for benchmark in BENCHMARKS:

        expected_n = (
            5303
            if benchmark == "panther_family"
            else 5323
        )

        x = audit_df[
            audit_df["benchmark"]
            == benchmark
        ]

        if not (
            x["n_genes"] == expected_n
        ).all():
            raise RuntimeError(
                f"{benchmark}: population "
                "size mismatch."
            )

    audit_path = (
        outroot
        / "split_audit.csv"
    )

    audit_df.to_csv(
        audit_path,
        index=False,
    )

    print("[5/5] Write protocol manifest")

    manifest = {
        "protocol":
            "MorphoShift_v5_paired",

        "performance_unblinded":
            False,

        "phenotype_columns_read":
            False,

        "canonical_gene_key":
            "Metadata_NCBI_Gene_ID",

        "paired_universe_n":
            5323,

        "panther_eligible_n":
            5303,

        "panther_missing_family_n":
            20,

        "directions": [
            "ORF_to_CRISPR",
            "CRISPR_to_ORF",
        ],

        "direction_split_policy":
            (
                "Both prediction directions use "
                "the identical train/val/test "
                "gene assignment."
            ),

        "ratios": RATIOS,

        "seeds": SEEDS,

        "benchmarks": {
            "random_gene": {
                "unit": "gene",
            },
            "homology30": {
                "unit": "cluster_id30",
            },
            "homology50": {
                "unit": "cluster_id50",
            },
            "panther_family": {
                "unit": "panther_family",
                "missing_family_policy":
                    "exclude",
            },
        },

        "paired_universe_source": {
            "path":
                str(pair_path.relative_to(root)),
            "sha256":
                sha256_file(pair_path),
        },

        "group_definition_sources":
            source_manifest,

        "split_algorithm": {
            "random_gene":
                (
                    "seeded gene-level random "
                    "permutation with largest-"
                    "remainder 80/10/10 counts"
                ),

            "grouped_benchmarks":
                (
                    "whole-group assignment; "
                    "256 deterministic randomized "
                    "greedy restarts minimizing "
                    "normalized deviation from "
                    "80/10/10 by gene count"
                ),

            "n_restarts":
                args.n_restarts,
        },

        "important":
            (
                "No morphology target X_* column "
                "was read while defining these "
                "splits."
            ),
    }

    with open(
        outroot / "split_protocol_manifest.json",
        "w",
    ) as f:
        json.dump(
            manifest,
            f,
            indent=2,
            ensure_ascii=False,
        )

    print()
    print("=== PAIRED SPLIT AUDIT ===")
    print(
        audit_df.to_string(index=False)
    )

    print()
    print("[DONE]", outroot)
    print(
        "[IMPORTANT] No morphology target "
        "column was read."
    )


if __name__ == "__main__":
    main()
