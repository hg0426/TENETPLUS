#!/usr/bin/env python3
"""Fast indirect-edge trimming for TENET SIF networks.

TENET-style trimming removes Source->Target edges when a stronger
Source->Middle->Target path exists. This implementation keeps the same input
and output contract, but replaces repeated list scans with dictionaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trim indirect edges from a TENET .sif network."
    )
    parser.add_argument("sif", type=Path, help="Input .sif file: Source TE Target.")
    parser.add_argument("threshold", type=float, help="Legacy indirect trim threshold.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path. Default: <input>.trimIndirect<threshold>.sif",
    )
    return parser.parse_args()


def default_output_path(sif_path: Path, threshold: float) -> Path:
    return sif_path.with_name(
        sif_path.name.replace(".sif", f".trimIndirect{threshold}.sif")
    )


def load_sif(sif_path: Path):
    source_order: list[str] = []
    targets_by_source: dict[str, list[str]] = {}
    tes_by_source: dict[str, list[float]] = {}

    with sif_path.open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 3:
                continue
            source, te_str, target = parts[0], parts[1], parts[2]
            if source not in targets_by_source:
                source_order.append(source)
                targets_by_source[source] = []
                tes_by_source[source] = []
            targets_by_source[source].append(target)
            tes_by_source[source].append(float(te_str))

    return source_order, targets_by_source, tes_by_source


def trim_indirect(
    source_order: list[str],
    targets_by_source: dict[str, list[str]],
    tes_by_source: dict[str, list[float]],
    threshold: float,
) -> dict[str, list[bool]]:
    source_set = set(source_order)
    target_te_by_source: dict[str, dict[str, float]] = {}
    for source in source_order:
        target_te: dict[str, float] = {}
        for target, te in zip(targets_by_source[source], tes_by_source[source]):
            if target not in target_te:
                target_te[target] = te
        target_te_by_source[source] = target_te
    indirect_by_source = {
        source: [False] * len(targets_by_source[source])
        for source in source_order
    }

    for source in source_order:
        source_targets = targets_by_source[source]
        source_tes = tes_by_source[source]
        source_indirect = indirect_by_source[source]

        for mid_idx, middle in enumerate(source_targets):
            if middle not in source_set:
                continue

            middle_target_te = target_te_by_source[middle]
            if not middle_target_te:
                continue

            source_middle_te = source_tes[mid_idx]
            for target_idx, target in enumerate(source_targets):
                if target_idx == mid_idx:
                    continue
                middle_target_te_value = middle_target_te.get(target)
                if middle_target_te_value is None:
                    continue
                cutoff = min(source_middle_te, middle_target_te_value) + threshold
                if source_tes[target_idx] < cutoff:
                    source_indirect[target_idx] = True

    return indirect_by_source


def write_trimmed(
    output_path: Path,
    source_order: list[str],
    targets_by_source: dict[str, list[str]],
    tes_by_source: dict[str, list[float]],
    indirect_by_source: dict[str, list[bool]],
) -> tuple[int, int]:
    kept = 0
    removed = 0
    with output_path.open("w") as handle:
        for source in source_order:
            targets = targets_by_source[source]
            tes = tes_by_source[source]
            indirect = indirect_by_source[source]
            for target, te, is_indirect in zip(targets, tes, indirect):
                if is_indirect:
                    removed += 1
                    continue
                kept += 1
                handle.write(f"{source}\t{te}\t{target}\n")
    return kept, removed


def main() -> None:
    args = parse_args()
    input_path = args.sif
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    output_path = args.output or default_output_path(input_path, args.threshold)
    start = perf_counter()
    source_order, targets_by_source, tes_by_source = load_sif(input_path)
    n_edges = sum(len(v) for v in targets_by_source.values())
    indirect_by_source = trim_indirect(
        source_order, targets_by_source, tes_by_source, args.threshold
    )
    kept, removed = write_trimmed(
        output_path, source_order, targets_by_source, tes_by_source, indirect_by_source
    )
    elapsed = perf_counter() - start
    print(
        "trim_indirect completed: "
        f"input={input_path} output={output_path} "
        f"sources={len(source_order)} edges={n_edges} kept={kept} "
        f"removed={removed} seconds={elapsed:.3f}"
    )


if __name__ == "__main__":
    main()
