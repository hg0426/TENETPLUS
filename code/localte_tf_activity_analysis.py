#!/usr/bin/env python3
"""Summarize timepoint LocalTE outdegree into group-specific TF activity.

The input is the compact output produced by ``make_GRN_localTE.py`` or
``trim_indirect_localte_times.py`` in ``outdegree`` mode:

    Time, Source, Outdegree

The script writes group-level TF activity tables and a group-enriched lineplot
for the top TFs per group. It is intentionally dataset-agnostic; groups can be
stages, cell types, clusters, conditions, lineages, or any other timepoint
annotation. Group labels can be supplied directly as a time-to-group table or
inferred from cell barcode prefixes in ``cell_gene_trsps.parquet``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

if "MPLCONFIGDIR" not in os.environ:
    _mpl_cache = tempfile.TemporaryDirectory(prefix="tenetplus_mpl_")
    os.environ["MPLCONFIGDIR"] = _mpl_cache.name

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


LINE_COLORS = [
    "#24364B",
    "#087783",
    "#C46B24",
    "#6B5FA7",
    "#2F8D46",
    "#B44765",
    "#6F6A2D",
    "#496FA8",
]

GROUP_COLORS = [
    "#2C7BB6",
    "#D7191C",
    "#FDAE61",
    "#1A9641",
    "#7B3294",
    "#A6611A",
    "#008C8C",
    "#666666",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdegree", required=True, type=Path, help="Parquet file or dataset with Time, Source, Outdegree columns.")
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for activity tables and plots.")
    parser.add_argument(
        "--time-group",
        "--time-stage",
        dest="time_group",
        type=Path,
        default=None,
        help=(
            "Optional TSV/CSV/parquet with Time or TimeInt plus a grouping column. "
            "--time-stage is kept as a backward-compatible alias."
        ),
    )
    parser.add_argument(
        "--group-column",
        default=None,
        help=(
            "Column in --time-group to use as the grouping label, such as Stage, "
            "celltype, cluster, condition, or lineage. If omitted, common names "
            "are auto-detected."
        ),
    )
    parser.add_argument(
        "--infer-group-from-cell-gene",
        "--infer-stage-from-cell-gene",
        dest="infer_group_from_cell_gene",
        type=Path,
        default=None,
        help=(
            "Infer group labels from cell columns in output/cell_gene_trsps.parquet. "
            "--infer-stage-from-cell-gene is kept as a backward-compatible alias."
        ),
    )
    parser.add_argument("--history-length", type=int, default=1, help="History length k used by TENETPLUS. Time 1 maps to cell column k.")
    parser.add_argument(
        "--group-regex",
        "--stage-regex",
        dest="group_regex",
        default=r"^([^_]+)",
        help="Regex used to extract a group from a cell barcode; group 1 is used when present.",
    )
    parser.add_argument(
        "--group-order",
        "--stage-order",
        dest="group_order",
        default="",
        help="Comma-separated group order. Default uses first appearance in pseudotime.",
    )
    parser.add_argument(
        "--group-label",
        default=None,
        help="Display name for the grouping dimension in plot titles. Default: group column name or Group.",
    )
    parser.add_argument("--max-time", type=int, default=None, help="Optional maximum TimeInt to include in the analysis.")
    parser.add_argument("--batch-size", type=int, default=1_000_000, help="Arrow scan batch size for outdegree aggregation.")
    parser.add_argument(
        "--top-per-group",
        "--top-per-stage",
        dest="top_per_group",
        type=int,
        default=3,
        help="Number of TFs selected per group for the lineplot.",
    )
    parser.add_argument("--rolling-window", type=int, default=501, help="Centered rolling window for plotted outdegree trajectories.")
    parser.add_argument("--candidate-top-n", type=int, default=0, help="Select group top TFs only from top N TFs by total outdegree; 0 means all.")
    parser.add_argument(
        "--group-specific-min-enrichment",
        "--stage-specific-min-enrichment",
        dest="group_specific_min_enrichment",
        type=float,
        default=1.25,
        help="Preferred minimum group mean / non-group mean.",
    )
    parser.add_argument(
        "--group-specific-min-fraction",
        "--stage-specific-min-fraction",
        dest="group_specific_min_fraction",
        type=float,
        default=0.15,
        help="Preferred minimum fraction of a TF's activity in the group.",
    )
    parser.add_argument(
        "--min-group-mean",
        "--min-stage-mean",
        dest="min_group_mean",
        type=float,
        default=0.0,
        help="Preferred minimum mean outdegree in the group.",
    )
    parser.add_argument("--plot-file", default="group_enriched_tf_outdegree_lineplots.png", help="Output plot filename.")
    parser.add_argument(
        "--group-background-alpha",
        "--stage-background-alpha",
        dest="group_background_alpha",
        type=float,
        default=0.24,
        help="Opacity for the group background bands in the lineplot.",
    )
    parser.add_argument("--dpi", type=int, default=220, help="Plot resolution.")
    return parser.parse_args()


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def _outdegree_dataset(path: Path) -> ds.Dataset:
    if not path.exists():
        raise FileNotFoundError(path)
    dataset = ds.dataset(str(path), format="parquet")
    required = {"Time", "Source", "Outdegree"}
    missing = sorted(required.difference(dataset.schema.names))
    if missing:
        raise ValueError(f"Outdegree dataset is missing required columns: {missing}")
    return dataset


def _iter_outdegree_batches(path: Path, columns: list[str], batch_size: int) -> Iterable[pd.DataFrame]:
    scanner = _outdegree_dataset(path).scanner(columns=columns, batch_size=max(1, int(batch_size)), use_threads=True)
    for batch in scanner.to_batches():
        if batch.num_rows:
            yield batch.to_pandas()


def _detect_group_column(df: pd.DataFrame, group_column: str | None) -> str:
    if group_column:
        if group_column not in df.columns:
            raise ValueError(f"--group-column {group_column!r} was not found in the time-group table")
        return group_column
    candidates = [
        "Group",
        "Stage",
        "stage",
        "CellType",
        "celltype",
        "cell_type",
        "Cell_Type",
        "cluster",
        "Cluster",
        "condition",
        "Condition",
        "lineage",
        "Lineage",
        "label",
        "Label",
    ]
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    non_time = [col for col in df.columns if col not in {"Time", "TimeInt"}]
    if len(non_time) == 1:
        return non_time[0]
    raise ValueError(
        "Could not auto-detect a grouping column. Provide --group-column "
        "with a column such as Stage, celltype, cluster, condition, or lineage."
    )


def _normalise_time_group(df: pd.DataFrame, group_column: str | None) -> pd.DataFrame:
    if "TimeInt" not in df.columns:
        if "Time" not in df.columns:
            raise ValueError("time-group table needs a TimeInt or Time column")
        df = df.copy()
        df["TimeInt"] = pd.to_numeric(df["Time"], errors="raise").astype(np.int64)
    group_col = _detect_group_column(df, group_column)
    out = df.copy()
    out["TimeInt"] = pd.to_numeric(out["TimeInt"], errors="raise").astype(np.int64)
    if bool(out[group_col].isna().any()) or bool(out[group_col].astype(str).str.strip().eq("").any()):
        raise ValueError(f"Grouping column {group_col!r} contains missing or empty labels.")
    out["Group"] = out[group_col].astype(str)
    out["GroupColumn"] = group_col
    conflicting = out.groupby("TimeInt", sort=False)["Group"].nunique()
    conflicting = conflicting.loc[conflicting > 1]
    if not conflicting.empty:
        examples = conflicting.index.astype(int).tolist()[:10]
        raise ValueError(f"TimeInt values map to multiple groups: {examples}")
    out = out.drop_duplicates("TimeInt").sort_values("TimeInt").reset_index(drop=True)
    leading_cols = ["TimeInt", "Group"]
    out = out.loc[:, leading_cols + [col for col in out.columns if col not in leading_cols]]
    if out.empty:
        raise ValueError("time-group table is empty")
    return out


def _scan_time_bounds(outdegree_path: Path, batch_size: int) -> tuple[int, int]:
    min_time: int | None = None
    max_time: int | None = None
    for batch in _iter_outdegree_batches(outdegree_path, ["Time"], batch_size):
        times = pd.to_numeric(batch["Time"], errors="raise").astype(np.int64)
        if times.empty:
            continue
        batch_min = int(times.min())
        batch_max = int(times.max())
        min_time = batch_min if min_time is None else min(min_time, batch_min)
        max_time = batch_max if max_time is None else max(max_time, batch_max)
    if min_time is None or max_time is None:
        raise ValueError(f"No Time values found in {outdegree_path}")
    return min_time, max_time


def _infer_time_group_from_cell_gene(
    cell_gene_path: Path,
    *,
    max_time: int,
    history_length: int,
    group_regex: str,
) -> pd.DataFrame:
    if not cell_gene_path.exists():
        raise FileNotFoundError(cell_gene_path)
    schema = ds.dataset(str(cell_gene_path), format="parquet").schema
    cell_cols = [name for name in schema.names if name != "__index_level_0__"]
    pattern = re.compile(group_regex)
    rows: list[dict[str, object]] = []
    for time_int in range(1, int(max_time) + 1):
        cell_idx = time_int + int(history_length) - 1
        if cell_idx >= len(cell_cols):
            raise ValueError(
                f"Cannot map TimeInt={time_int} to cell column {cell_idx}; "
                f"only {len(cell_cols)} cell columns are available."
            )
        cell = str(cell_cols[cell_idx])
        match = pattern.search(cell)
        if match:
            group = match.group(1) if match.groups() else match.group(0)
        else:
            group = cell
        rows.append({"TimeInt": time_int, "Cell": cell, "Group": str(group)})
    return pd.DataFrame(rows)


def load_time_group(args: argparse.Namespace) -> pd.DataFrame:
    if args.time_group and args.infer_group_from_cell_gene:
        raise ValueError("Use either --time-group or --infer-group-from-cell-gene, not both.")
    if args.time_group:
        group = _normalise_time_group(_read_table(args.time_group), args.group_column)
    elif args.infer_group_from_cell_gene:
        _, scanned_max_time = _scan_time_bounds(args.outdegree, args.batch_size)
        max_time = int(args.max_time) if args.max_time is not None else scanned_max_time
        group = _infer_time_group_from_cell_gene(
            args.infer_group_from_cell_gene,
            max_time=max_time,
            history_length=max(1, int(args.history_length)),
            group_regex=args.group_regex,
        )
        group = _normalise_time_group(group, "Group")
    else:
        raise ValueError("Provide --time-group or --infer-group-from-cell-gene.")

    if args.max_time is not None:
        group = group.loc[group["TimeInt"] <= int(args.max_time)].copy()
    if group.empty:
        raise ValueError("No group labels remain after applying --max-time.")
    return group


def _ordered_groups(group_map: pd.DataFrame, group_order_arg: str) -> list[str]:
    if group_order_arg.strip():
        requested = [item.strip() for item in group_order_arg.split(",") if item.strip()]
        observed = set(group_map["Group"].astype(str))
        missing = [group for group in requested if group not in observed]
        if missing:
            raise ValueError(f"--group-order contains groups not found in the time-group map: {missing}")
        extras = [group for group in group_map["Group"].astype(str).drop_duplicates().tolist() if group not in requested]
        return requested + extras
    return group_map["Group"].astype(str).drop_duplicates().tolist()


def _combine_max(existing: pd.Series | None, piece: pd.Series) -> pd.Series:
    if existing is None:
        return piece
    both = pd.concat([existing.rename("old"), piece.rename("new")], axis=1)
    return both.max(axis=1)


def aggregate_outdegree(
    *,
    outdegree_path: Path,
    group_map: pd.DataFrame,
    batch_size: int,
    max_time: int | None = None,
) -> tuple[pd.DataFrame, pd.Series, dict[str, dict[str, object]]]:
    group_by_time = group_map.set_index("TimeInt")["Group"].astype(str)
    group_sum: pd.Series | None = None
    group_count: pd.Series | None = None
    group_max: pd.Series | None = None
    source_total: pd.Series | None = None
    source_peak: dict[str, dict[str, object]] = {}
    rows_seen = 0

    for batch in _iter_outdegree_batches(outdegree_path, ["Time", "Source", "Outdegree"], batch_size):
        work = batch.copy()
        work["TimeInt"] = pd.to_numeric(work["Time"], errors="raise").astype(np.int64)
        if max_time is not None:
            work = work.loc[work["TimeInt"] <= int(max_time)].copy()
        if work.empty:
            continue
        work["Group"] = work["TimeInt"].map(group_by_time)
        missing = work["Group"].isna()
        if bool(missing.any()):
            examples = sorted(work.loc[missing, "TimeInt"].drop_duplicates().head(10).astype(int).tolist())
            raise ValueError(f"Group map does not cover TimeInt values such as {examples}")
        work["Group"] = work["Group"].astype(str)
        work["Source"] = work["Source"].astype(str)
        work["Outdegree"] = pd.to_numeric(work["Outdegree"], errors="raise").astype(np.float64)
        rows_seen += len(work)

        total_piece = work.groupby("Source", sort=False)["Outdegree"].sum()
        source_total = total_piece if source_total is None else source_total.add(total_piece, fill_value=0.0)

        grouped = work.groupby(["Source", "Group"], sort=False)["Outdegree"]
        sum_piece = grouped.sum()
        count_piece = grouped.size().astype(np.int64)
        max_piece = grouped.max()
        group_sum = sum_piece if group_sum is None else group_sum.add(sum_piece, fill_value=0.0)
        group_count = count_piece if group_count is None else group_count.add(count_piece, fill_value=0)
        group_max = _combine_max(group_max, max_piece)

        idx = work.groupby("Source", sort=False)["Outdegree"].idxmax()
        for row in work.loc[idx, ["Source", "TimeInt", "Group", "Outdegree"]].itertuples(index=False):
            source = str(row.Source)
            value = float(row.Outdegree)
            current = source_peak.get(source)
            if current is None or value > float(current["peak_outdegree"]):
                source_peak[source] = {
                    "peak_time": int(row.TimeInt),
                    "peak_group": str(row.Group),
                    "peak_outdegree": value,
                }

    if rows_seen == 0 or source_total is None or group_sum is None or group_count is None:
        raise ValueError(f"No rows found in outdegree dataset: {outdegree_path}")

    if group_max is None:
        group_max = pd.Series(dtype=float)

    group_counts = group_map.groupby("Group", sort=False)["TimeInt"].nunique().astype(np.int64)
    total_timepoints = int(group_counts.sum())
    sources = source_total.sort_values(ascending=False).index.astype(str).tolist()
    groups = group_counts.index.astype(str).tolist()
    full_index = pd.MultiIndex.from_product([sources, groups], names=["Source", "Group"])

    stats = pd.DataFrame(index=full_index).reset_index()
    stats["group_outdegree"] = group_sum.reindex(full_index, fill_value=0.0).to_numpy(dtype=float)
    stats["active_timepoints"] = group_count.reindex(full_index, fill_value=0).to_numpy(dtype=np.int64)
    stats["max_outdegree"] = group_max.reindex(full_index, fill_value=0.0).to_numpy(dtype=float)
    stats["group_timepoints"] = stats["Group"].map(group_counts).astype(np.int64)
    stats["tf_total_outdegree"] = stats["Source"].map(source_total).astype(float)
    stats["group_mean_outdegree"] = stats["group_outdegree"] / stats["group_timepoints"].clip(lower=1)
    stats["group_fraction_of_tf"] = np.divide(
        stats["group_outdegree"],
        stats["tf_total_outdegree"],
        out=np.zeros(len(stats), dtype=float),
        where=stats["tf_total_outdegree"].to_numpy(dtype=float) > 0,
    )

    non_group_timepoints = total_timepoints - stats["group_timepoints"]
    non_group_outdegree = stats["tf_total_outdegree"] - stats["group_outdegree"]
    stats["non_group_mean_outdegree"] = np.divide(
        non_group_outdegree,
        non_group_timepoints.clip(lower=1),
        out=np.zeros(len(stats), dtype=float),
        where=non_group_timepoints.to_numpy(dtype=float) > 0,
    )
    global_mean = stats["tf_total_outdegree"] / max(1, total_timepoints)
    eps = 1e-9
    stats["group_enrichment_vs_non_group"] = (stats["group_mean_outdegree"] + eps) / (
        stats["non_group_mean_outdegree"] + eps
    )
    stats["group_enrichment_vs_global"] = (stats["group_mean_outdegree"] + eps) / (global_mean + eps)

    peak_df = pd.DataFrame.from_dict(source_peak, orient="index")
    peak_df.index.name = "Source"
    peak_df = peak_df.reset_index()
    stats = stats.merge(peak_df, on="Source", how="left")
    stats["rank_within_group"] = (
        stats.sort_values(
            ["Group", "group_enrichment_vs_non_group", "group_fraction_of_tf", "group_mean_outdegree"],
            ascending=[True, False, False, False],
        )
        .groupby("Group")
        .cumcount()
        + 1
    )
    return stats, source_total.sort_values(ascending=False), source_peak


def select_group_tfs(
    stats: pd.DataFrame,
    *,
    group_order: list[str],
    source_total: pd.Series,
    top_per_group: int,
    candidate_top_n: int,
    min_enrichment: float,
    min_fraction: float,
    min_group_mean: float,
) -> pd.DataFrame:
    if candidate_top_n > 0:
        candidates = set(source_total.head(candidate_top_n).index.astype(str))
        pool = stats.loc[stats["Source"].isin(candidates)].copy()
    else:
        pool = stats.copy()

    sort_cols = [
        "group_enrichment_vs_non_group",
        "group_fraction_of_tf",
        "group_mean_outdegree",
        "group_outdegree",
        "tf_total_outdegree",
    ]
    rows: list[pd.DataFrame] = []
    for group in group_order:
        part = pool.loc[pool["Group"] == group].copy()
        if part.empty:
            continue
        part = part.sort_values(sort_cols, ascending=[False, False, False, False, False])
        passed = part.loc[
            (part["group_enrichment_vs_non_group"] >= min_enrichment)
            & (part["group_fraction_of_tf"] >= min_fraction)
            & (part["group_mean_outdegree"] >= min_group_mean)
        ].copy()
        passed["passed_specificity_filter"] = True
        if len(passed) < top_per_group:
            fill = part.loc[~part["Source"].isin(passed["Source"])].head(top_per_group - len(passed)).copy()
            fill["passed_specificity_filter"] = False
            chosen = pd.concat([passed, fill], ignore_index=True)
        else:
            chosen = passed.head(top_per_group)
        chosen = chosen.head(top_per_group).copy()
        chosen["group_plot_rank"] = np.arange(1, len(chosen) + 1)
        rows.append(chosen)
    if not rows:
        raise ValueError("No TFs were selected for plotting.")
    return pd.concat(rows, ignore_index=True)


def build_selected_timecourse(
    *,
    outdegree_path: Path,
    group_map: pd.DataFrame,
    selected_sources: list[str],
    batch_size: int,
    rolling_window: int,
) -> pd.DataFrame:
    time_index = pd.Index(group_map["TimeInt"].astype(np.int64), name="TimeInt")
    selected_sources = list(dict.fromkeys(str(source) for source in selected_sources))
    wide = pd.DataFrame(0.0, index=time_index, columns=selected_sources, dtype=np.float32)
    selected_set = set(selected_sources)

    for batch in _iter_outdegree_batches(outdegree_path, ["Time", "Source", "Outdegree"], batch_size):
        work = batch.copy()
        work["Source"] = work["Source"].astype(str)
        work = work.loc[work["Source"].isin(selected_set)]
        if work.empty:
            continue
        work["TimeInt"] = pd.to_numeric(work["Time"], errors="raise").astype(np.int64)
        work["Outdegree"] = pd.to_numeric(work["Outdegree"], errors="raise").astype(np.float32)
        grouped = work.groupby(["TimeInt", "Source"], sort=False)["Outdegree"].sum().reset_index()
        pivot = grouped.pivot(index="TimeInt", columns="Source", values="Outdegree")
        pivot = pivot.reindex(index=wide.index, columns=wide.columns, fill_value=0.0).fillna(0.0).astype(np.float32)
        wide += pivot

    min_periods = max(1, min(int(rolling_window), len(wide)) // 10)
    rolled = (
        wide.rolling(window=max(1, int(rolling_window)), center=True, min_periods=min_periods)
        .mean()
        .bfill()
        .ffill()
        .fillna(0.0)
        .astype(np.float32)
    )

    raw_long = wide.reset_index().melt(id_vars="TimeInt", var_name="Source", value_name="Outdegree")
    roll_long = rolled.reset_index().melt(id_vars="TimeInt", var_name="Source", value_name="SmoothedOutdegree")
    long = raw_long.merge(roll_long, on=["TimeInt", "Source"], how="inner")
    long = long.merge(group_map[["TimeInt", "Group"]], on="TimeInt", how="left")
    return long.loc[:, ["TimeInt", "Group", "Source", "Outdegree", "SmoothedOutdegree"]]


def _group_spans(group_map: pd.DataFrame) -> list[tuple[int, int, str]]:
    ordered = group_map.sort_values("TimeInt")
    spans: list[tuple[int, int, str]] = []
    start = int(ordered["TimeInt"].iloc[0])
    prev_time = start
    prev_group = str(ordered["Group"].iloc[0])
    for row in ordered.iloc[1:].itertuples(index=False):
        time_int = int(row.TimeInt)
        group = str(row.Group)
        if group != prev_group:
            spans.append((start, prev_time, prev_group))
            start = time_int
            prev_group = group
        prev_time = time_int
    spans.append((start, prev_time, prev_group))
    return spans


def plot_group_lineplots(
    *,
    timecourse: pd.DataFrame,
    selected: pd.DataFrame,
    group_map: pd.DataFrame,
    group_order: list[str],
    outpath: Path,
    top_per_group: int,
    rolling_window: int,
    dpi: int,
    group_label: str,
    group_background_alpha: float,
) -> None:
    groups_with_rows = [group for group in group_order if group in set(selected["Group"].astype(str))]
    if not groups_with_rows:
        raise ValueError("No selected groups are available for plotting.")

    group_color_map = {group: GROUP_COLORS[i % len(GROUP_COLORS)] for i, group in enumerate(group_order)}
    spans = _group_spans(group_map)
    span_groups = list(dict.fromkeys(span_group for _, _, span_group in spans))
    group_handles = [
        Patch(
            facecolor=group_color_map.get(group, "#CCCCCC"),
            edgecolor="#333333",
            linewidth=0.6,
            alpha=0.85,
            label=group,
        )
        for group in group_order
        if group in span_groups
    ]
    nrows = len(groups_with_rows)
    fig_height = max(2.3 * nrows, 4.0)
    fig, axes = plt.subplots(nrows, 1, figsize=(14, fig_height), sharex=True, squeeze=False)
    axes_flat = axes[:, 0]

    for ax, group in zip(axes_flat, groups_with_rows):
        for start, end, span_group in spans:
            ax.axvspan(
                start - 0.5,
                end + 0.5,
                color=group_color_map.get(span_group, "#EEEEEE"),
                alpha=group_background_alpha,
                linewidth=0,
            )
        group_selected = selected.loc[selected["Group"] == group].sort_values("group_plot_rank")
        for color_idx, source in enumerate(group_selected["Source"].astype(str)):
            part = timecourse.loc[timecourse["Source"].astype(str) == source].sort_values("TimeInt")
            ax.plot(
                part["TimeInt"].to_numpy(dtype=int),
                part["SmoothedOutdegree"].to_numpy(dtype=float),
                label=source,
                color=LINE_COLORS[color_idx % len(LINE_COLORS)],
                linewidth=1.5,
            )
        ax.set_ylabel(group)
        ax.grid(True, axis="y", color="#E5E7EB", linewidth=0.7)
        ax.grid(False, axis="x")
        ax.legend(loc="upper right", ncol=min(top_per_group, 4), frameon=False, fontsize=9)
        ax.tick_params(axis="both", labelsize=9)

    axes_flat[-1].set_xlabel("Pseudotime index")
    fig.supylabel("LocalTE outdegree", x=0.015, fontsize=11)
    title = f"{group_label}-enriched TF LocalTE outdegree lineplots (top {top_per_group} per group, rolling window={rolling_window})"
    fig.suptitle(title, fontsize=15, y=0.995)
    if group_handles:
        fig.legend(
            handles=group_handles,
            title=f"{group_label} bands",
            loc="upper center",
            bbox_to_anchor=(0.5, 0.965),
            ncol=min(len(group_handles), 8),
            frameon=False,
            fontsize=9.5,
            title_fontsize=10,
            handlelength=1.8,
            handleheight=0.9,
            columnspacing=1.5,
        )
    fig.tight_layout(rect=[0.03, 0.025, 0.995, 0.925])
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    group_map = load_time_group(args)
    group_order = _ordered_groups(group_map, args.group_order)
    group_label = args.group_label or str(group_map["GroupColumn"].iloc[0] if "GroupColumn" in group_map.columns else "Group")
    group_map.to_csv(args.outdir / "time_to_group_map.tsv", sep="\t", index=False)

    stats, source_total, _ = aggregate_outdegree(
        outdegree_path=args.outdegree,
        group_map=group_map,
        batch_size=args.batch_size,
        max_time=args.max_time,
    )
    stats.to_csv(args.outdir / "tf_group_activity.tsv", sep="\t", index=False)
    source_total.rename("TotalOutdegree").to_csv(args.outdir / "top_by_total_outdegree.tsv", sep="\t")

    selected = select_group_tfs(
        stats,
        group_order=group_order,
        source_total=source_total,
        top_per_group=max(1, int(args.top_per_group)),
        candidate_top_n=max(0, int(args.candidate_top_n)),
        min_enrichment=float(args.group_specific_min_enrichment),
        min_fraction=float(args.group_specific_min_fraction),
        min_group_mean=float(args.min_group_mean),
    )
    selected.to_csv(args.outdir / "group_specific_top_tfs.tsv", sep="\t", index=False)

    timecourse = build_selected_timecourse(
        outdegree_path=args.outdegree,
        group_map=group_map,
        selected_sources=selected["Source"].astype(str).tolist(),
        batch_size=args.batch_size,
        rolling_window=max(1, int(args.rolling_window)),
    )
    timecourse.to_parquet(args.outdir / "selected_tf_timecourse_outdegree_smoothed.parquet", index=False)

    plot_path = args.outdir / args.plot_file
    plot_group_lineplots(
        timecourse=timecourse,
        selected=selected,
        group_map=group_map,
        group_order=group_order,
        outpath=plot_path,
        top_per_group=max(1, int(args.top_per_group)),
        rolling_window=max(1, int(args.rolling_window)),
        dpi=max(72, int(args.dpi)),
        group_label=group_label,
        group_background_alpha=min(max(float(args.group_background_alpha), 0.0), 1.0),
    )

    summary = {
        "outdegree": str(args.outdegree),
        "outdir": str(args.outdir),
        "group_label": group_label,
        "n_timepoints": int(group_map["TimeInt"].nunique()),
        "n_groups": int(group_map["Group"].nunique()),
        "n_sources": int(source_total.shape[0]),
        "selected_sources": sorted(set(selected["Source"].astype(str))),
        "top_per_group": int(args.top_per_group),
        "rolling_window": int(args.rolling_window),
        "group_background_alpha": min(max(float(args.group_background_alpha), 0.0), 1.0),
        "plot": str(plot_path),
        "specificity_filters": {
            "group_specific_min_enrichment": float(args.group_specific_min_enrichment),
            "group_specific_min_fraction": float(args.group_specific_min_fraction),
            "min_group_mean": float(args.min_group_mean),
        },
    }
    (args.outdir / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
