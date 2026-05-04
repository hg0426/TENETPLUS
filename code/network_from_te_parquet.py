#!/usr/bin/env python3
"""Build sparse network files from a large TENET TE parquet table.

The input is expected to have Source, Target, and TE columns.  The script uses
DuckDB so a full all-feature TE table does not need to be loaded into pandas.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import norm


def _sql_path(path: Path) -> str:
    return str(path.resolve()).replace("'", "''")


def _write_sif(con: duckdb.DuckDBPyConnection, parquet_path: Path, sif_path: Path) -> None:
    con.execute(
        f"""
        COPY (
            SELECT
                CAST(Source AS VARCHAR) AS Source,
                CAST(TE AS VARCHAR) AS TE,
                CAST(Target AS VARCHAR) AS Target
            FROM read_parquet('{_sql_path(parquet_path)}')
        )
        TO '{_sql_path(sif_path)}' (FORMAT CSV, DELIMITER '\t', HEADER false)
        """
    )


def _fdr_threshold(
    con: duckdb.DuckDBPyConnection,
    input_path: Path,
    alpha: float,
    te_cutoff: float,
    batch_size: int = 1_000_000,
) -> dict[str, float | int]:
    """Find the exact global z-score + BH-FDR TE cutoff.

    TENET's GRN step first filters TE > te_cutoff, computes a global z-score
    over those retained TE values, converts high-z edges to one-sided p-values,
    then applies BH-FDR. Since p-values are monotone in TE, we stream TE values
    sorted descending and test the exact BH rank boundary without loading the
    full parquet into pandas.
    """
    m, mean, std, min_te, max_te = con.execute(
        f"""
        SELECT COUNT(*), AVG(TE), STDDEV_POP(TE), MIN(TE), MAX(TE)
        FROM read_parquet('{_sql_path(input_path)}')
        WHERE TE > {float(te_cutoff)}
        """
    ).fetchone()
    m = int(m or 0)
    if m < 2 or std is None or float(std) <= 0:
        return {
            "n_tested": m,
            "mean": float(mean or 0.0),
            "std": float(std or 0.0),
            "threshold": float("inf"),
            "n_edges": 0,
            "p_value": 1.0,
            "bh_bound": 0.0,
            "rank": 0,
        }

    mean_f = float(mean)
    std_f = float(std)
    best_rank = 0
    best_threshold = None
    best_p_value = 1.0
    best_bh_bound = 0.0
    offset = 0
    reader = con.execute(
        f"""
        SELECT TE
        FROM read_parquet('{_sql_path(input_path)}')
        WHERE TE > {float(te_cutoff)}
        ORDER BY TE DESC
        """
    ).fetch_record_batch(rows_per_batch=max(1, int(batch_size)))
    for batch in reader:
        te_values = batch.column(0).to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        if te_values.size == 0:
            continue
        ranks = np.arange(offset + 1, offset + te_values.size + 1, dtype=np.float64)
        p_values = norm.sf((te_values - mean_f) / std_f)
        bh_bounds = float(alpha) * ranks / float(m)
        ok = p_values <= bh_bounds
        if np.any(ok):
            last_ok = int(np.flatnonzero(ok)[-1])
            best_rank = int(offset + last_ok + 1)
            best_threshold = float(te_values[last_ok])
            best_p_value = float(p_values[last_ok])
            best_bh_bound = float(bh_bounds[last_ok])
        offset += int(te_values.size)

    if best_threshold is None:
        return {
            "n_tested": m,
            "mean": float(mean),
            "std": float(std),
            "threshold": float("inf"),
            "n_edges": 0,
            "p_value": 1.0,
            "bh_bound": 0.0,
            "rank": 0,
        }
    # Include all ties at the selected TE threshold, matching BH's p-value
    # threshold semantics rather than an arbitrary row_number cutoff.
    n_edges = con.execute(
        f"""
        SELECT COUNT(*)
        FROM read_parquet('{_sql_path(input_path)}')
        WHERE TE > {float(te_cutoff)} AND TE >= {best_threshold}
        """
    ).fetchone()[0]
    return {
        "n_tested": m,
        "mean": float(mean),
        "std": float(std),
        "threshold": float(best_threshold),
        "n_edges": int(n_edges),
        "p_value": float(best_p_value),
        "bh_bound": float(best_bh_bound),
        "rank": int(best_rank),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create network edge tables/SIF files from TENET TE parquet output."
    )
    parser.add_argument("input", type=Path, help="Input TE parquet with Source, Target, TE columns.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for output files. Default: input parquet stem, e.g. TE_TF_GN.",
    )
    parser.add_argument("--topk-per-target", type=int, default=0, help="Keep top K incoming sources per target.")
    parser.add_argument(
        "--fdr",
        type=float,
        dest="fdr",
        default=None,
        help="Apply global z-score + BH-FDR threshold and export <prefix>_fdr*. files.",
    )
    parser.add_argument(
        "--te-cutoff",
        type=float,
        dest="te_cutoff",
        default=0.0,
        help="TE cutoff used before global z-score/FDR.",
    )
    parser.add_argument(
        "--quantile",
        type=float,
        default=None,
        help="Global TE quantile threshold to export, e.g. 0.9999 keeps the top 0.01%%.",
    )
    parser.add_argument("--min-te", type=float, default=None, help="Optional absolute TE cutoff.")
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1), help="DuckDB threads.")
    parser.add_argument("--memory-limit", default=None, help="Optional DuckDB memory limit, e.g. 32GB.")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={int(args.threads)}")
    if args.memory_limit:
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")

    input_path = args.input.resolve()
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    output_prefix = args.output_prefix or input_path.stem

    summary_path = args.outdir / "network_summary.txt"
    with summary_path.open("w") as fh:
        fh.write(f"input\t{input_path}\n")
        fh.write(f"threads\t{args.threads}\n")
        if args.memory_limit:
            fh.write(f"memory_limit\t{args.memory_limit}\n")

        stats = con.execute(
            f"""
            SELECT
                COUNT(*) AS n_edges_raw,
                COUNT(DISTINCT Source) AS n_sources,
                COUNT(DISTINCT Target) AS n_targets,
                MIN(TE) AS min_te,
                MAX(TE) AS max_te,
                AVG(TE) AS mean_te
            FROM read_parquet('{_sql_path(input_path)}')
            """
        ).fetchone()
        fh.write(
            "raw_stats\t"
            f"n_edges={stats[0]}\t"
            f"n_sources={stats[1]}\t"
            f"n_targets={stats[2]}\t"
            f"min_te={stats[3]}\t"
            f"max_te={stats[4]}\t"
            f"mean_te={stats[5]}\n"
        )

        if args.quantile is not None:
            q = float(args.quantile)
            q_value = con.execute(
                f"SELECT quantile_cont(TE, {q}) FROM read_parquet('{_sql_path(input_path)}')"
            ).fetchone()[0]
            quantile_edges = args.outdir / f"{output_prefix}_global_q{q:g}.parquet"
            con.execute(
                f"""
                COPY (
                    SELECT
                        CAST(Source AS VARCHAR) AS Source,
                        CAST(Target AS VARCHAR) AS Target,
                        TE
                    FROM read_parquet('{_sql_path(input_path)}')
                    WHERE TE >= {float(q_value)}
                    ORDER BY TE DESC
                )
                TO '{_sql_path(quantile_edges)}' (FORMAT PARQUET)
                """
            )
            quantile_sif = args.outdir / f"{output_prefix}_global_q{q:g}.sif"
            _write_sif(con, quantile_edges, quantile_sif)
            n_q = con.execute(f"SELECT COUNT(*) FROM read_parquet('{_sql_path(quantile_edges)}')").fetchone()[0]
            fh.write(
                f"global_quantile\tq={q}\tthreshold={q_value}\t"
                f"n_edges={n_q}\tparquet={quantile_edges.name}\tsif={quantile_sif.name}\n"
            )

        if args.fdr is not None:
            alpha = float(args.fdr)
            cutoff = float(args.te_cutoff)
            fdr_stats = _fdr_threshold(con, input_path, alpha=alpha, te_cutoff=cutoff)
            fdr_edges = args.outdir / f"{output_prefix}_fdr{alpha:g}.parquet"
            threshold = float(fdr_stats["threshold"])
            if int(fdr_stats["n_edges"]) > 0:
                con.execute(
                    f"""
                    COPY (
                        SELECT
                            CAST(Source AS VARCHAR) AS Source,
                            CAST(Target AS VARCHAR) AS Target,
                            TE
                        FROM read_parquet('{_sql_path(input_path)}')
                        WHERE TE > {cutoff} AND TE >= {threshold}
                        ORDER BY TE DESC
                    )
                    TO '{_sql_path(fdr_edges)}' (FORMAT PARQUET)
                    """
                )
                fdr_sif = args.outdir / f"{output_prefix}_fdr{alpha:g}.sif"
                _write_sif(con, fdr_edges, fdr_sif)
            else:
                fdr_sif = args.outdir / f"{output_prefix}_fdr{alpha:g}.sif"
                con.execute(
                    f"""
                    COPY (
                        SELECT
                            CAST(NULL AS VARCHAR) AS Source,
                            CAST(NULL AS VARCHAR) AS Target,
                            CAST(NULL AS DOUBLE) AS TE
                        WHERE false
                    )
                    TO '{_sql_path(fdr_edges)}' (FORMAT PARQUET)
                    """
                )
                fdr_sif.write_text("")
            fh.write(
                f"fdr\talpha={alpha}\tte_cutoff={cutoff}\t"
                f"n_tested={fdr_stats['n_tested']}\tmean={fdr_stats['mean']}\tstd={fdr_stats['std']}\t"
                f"threshold={fdr_stats['threshold']}\tn_edges={fdr_stats['n_edges']}\t"
                f"bh_rank={fdr_stats.get('rank', 0)}\t"
                f"p_value_at_threshold={fdr_stats['p_value']}\tbh_bound={fdr_stats['bh_bound']}\t"
                f"parquet={fdr_edges.name}\tsif={fdr_sif.name}\n"
            )

        if args.min_te is not None:
            min_te = float(args.min_te)
            min_edges = args.outdir / f"{output_prefix}_minTE_{min_te:g}.parquet"
            con.execute(
                f"""
                COPY (
                    SELECT
                        CAST(Source AS VARCHAR) AS Source,
                        CAST(Target AS VARCHAR) AS Target,
                        TE
                    FROM read_parquet('{_sql_path(input_path)}')
                    WHERE TE >= {min_te}
                    ORDER BY TE DESC
                )
                TO '{_sql_path(min_edges)}' (FORMAT PARQUET)
                """
            )
            min_sif = args.outdir / f"{output_prefix}_minTE_{min_te:g}.sif"
            _write_sif(con, min_edges, min_sif)
            n_min = con.execute(f"SELECT COUNT(*) FROM read_parquet('{_sql_path(min_edges)}')").fetchone()[0]
            fh.write(
                f"min_te\tthreshold={min_te}\tn_edges={n_min}\t"
                f"parquet={min_edges.name}\tsif={min_sif.name}\n"
            )

        if int(args.topk_per_target) > 0:
            k = int(args.topk_per_target)
            topk_edges = args.outdir / f"{output_prefix}_top{k}_per_target.parquet"
            con.execute(
                f"""
                COPY (
                    SELECT Source, Target, TE, rank_in_target
                    FROM (
                        SELECT
                            CAST(Source AS VARCHAR) AS Source,
                            CAST(Target AS VARCHAR) AS Target,
                            TE,
                            ROW_NUMBER() OVER (
                                PARTITION BY Target
                                ORDER BY TE DESC, Source
                            ) AS rank_in_target
                        FROM read_parquet('{_sql_path(input_path)}')
                    )
                    WHERE rank_in_target <= {k}
                    ORDER BY Target, rank_in_target
                )
                TO '{_sql_path(topk_edges)}' (FORMAT PARQUET)
                """
            )
            topk_sif = args.outdir / f"{output_prefix}_top{k}_per_target.sif"
            _write_sif(con, topk_edges, topk_sif)
            n_topk = con.execute(f"SELECT COUNT(*) FROM read_parquet('{_sql_path(topk_edges)}')").fetchone()[0]
            fh.write(
                f"topk_per_target\tk={k}\tn_edges={n_topk}\t"
                f"parquet={topk_edges.name}\tsif={topk_sif.name}\n"
            )

    print(f"Wrote network outputs to {args.outdir}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
