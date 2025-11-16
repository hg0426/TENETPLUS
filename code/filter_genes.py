from __future__ import annotations

import argparse
from typing import List

import pandas as pd


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter genes from an expression matrix by name prefix."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input Parquet matrix path (cells x features).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output Parquet matrix path with filtered columns.",
    )
    parser.add_argument(
        "--exclude_prefixes",
        default="",
        help="Comma-separated list of prefixes to drop (e.g. 'RPS,RPL,MRPS,MRPL,MT-').",
    )
    return parser.parse_args(argv)


def filter_by_prefix(
    input_path: str,
    output_path: str,
    exclude_prefixes: List[str],
) -> None:
    df = pd.read_parquet(input_path)
    prefixes = [p for p in (s.strip() for s in exclude_prefixes) if p]
    if not prefixes:
        df.to_parquet(output_path)
        return
    keep_cols = [
        col
        for col in df.columns
        if not any(str(col).startswith(p) for p in prefixes)
    ]
    if not keep_cols:
        raise SystemExit(
            "No columns remain after applying gene filters; "
            "adjust TENET_GENE_FILTER or prefixes."
        )
    df.loc[:, keep_cols].to_parquet(output_path)
    removed = len(df.columns) - len(keep_cols)
    print(
        f"[filter_genes] input columns={len(df.columns)} "
        f"removed={removed} remaining={len(keep_cols)}"
    )


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)
    prefixes = args.exclude_prefixes.split(",") if args.exclude_prefixes else []
    filter_by_prefix(args.input, args.output, prefixes)


if __name__ == "__main__":
    main()

