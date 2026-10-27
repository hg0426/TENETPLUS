#!/usr/bin/env python3
"""Query a chunked LocalTE timepoint GRN parquet dataset."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import duckdb

from code.path_utils import coerce_input_path, coerce_output_path


def _sql_list(values: Optional[str]) -> Optional[list[str]]:
    if values is None or str(values).strip() == "":
        return None
    return [item.strip() for item in str(values).split(",") if item.strip()]


def _time_predicate(time_expr: str, time_filter: Optional[str]) -> tuple[str, list[object]]:
    if not time_filter:
        return "", []
    predicates: list[str] = []
    params: list[object] = []
    for part in str(time_filter).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            left, right = part.split("-", 1)
            predicates.append(f"try_cast({time_expr} as integer) between ? and ?")
            params.extend([int(left), int(right)])
        else:
            predicates.append(f"try_cast({time_expr} as integer) = ?")
            params.append(int(part))
    if not predicates:
        return "", []
    return "(" + " or ".join(predicates) + ")", params


def _sql_string_literal(value: Path) -> str:
    """DuckDB COPY target paths cannot be parameter-bound in all versions."""
    return "'" + str(value).replace("'", "''") + "'"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", required=True, help="Dataset directory or glob, e.g. localTE_GRN.dataset/*.parquet")
    p.add_argument("--output", help="Output parquet/csv/tsv. Omit for count/preview only.")
    p.add_argument("--source", help="Comma-separated source names")
    p.add_argument("--target", help="Comma-separated target names")
    p.add_argument("--time", help="Comma-separated timepoints/ranges, e.g. 1,5,100-200")
    p.add_argument("--columns", default="Source,TE,Target,Time", help="Comma-separated columns to select")
    p.add_argument("--limit", type=int, help="Optional row limit")
    p.add_argument("--count-only", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    dataset_raw = Path(args.dataset)
    if any(ch in str(args.dataset) for ch in "*?[]"):
        dataset_glob = str(dataset_raw)
    else:
        dataset_path = coerce_input_path(args.dataset)
        dataset_glob = str(dataset_path / "*.parquet") if dataset_path.is_dir() else str(dataset_path)

    columns = [col.strip() for col in args.columns.split(",") if col.strip()]
    if not columns:
        raise ValueError("--columns must select at least one column")
    select_sql = ", ".join(columns)
    where: list[str] = []
    params: list[object] = []
    sources = _sql_list(args.source)
    targets = _sql_list(args.target)
    if sources:
        where.append("Source in (" + ",".join(["?"] * len(sources)) + ")")
        params.extend(sources)
    if targets:
        where.append("Target in (" + ",".join(["?"] * len(targets)) + ")")
        params.extend(targets)
    time_clause, time_params = _time_predicate("Time", args.time)
    if time_clause:
        where.append(time_clause)
        params.extend(time_params)
    where_sql = (" where " + " and ".join(where)) if where else ""
    limit_sql = f" limit {int(args.limit)}" if args.limit else ""

    con = duckdb.connect()
    if args.count_only:
        sql = f"select count(*) as n_rows from read_parquet(?){where_sql}"
        print(con.execute(sql, [dataset_glob, *params]).df().to_string(index=False))
        return

    sql = f"select {select_sql} from read_parquet(?){where_sql}{limit_sql}"
    if args.output:
        output = coerce_output_path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output_sql = _sql_string_literal(output)
        if output.suffix.lower() == ".csv":
            copy_sql = f"copy ({sql}) to {output_sql} (header, delimiter ',')"
        elif output.suffix.lower() in {".tsv", ".txt"}:
            copy_sql = f"copy ({sql}) to {output_sql} (header, delimiter '\t')"
        else:
            copy_sql = f"copy ({sql}) to {output_sql} (format parquet, compression snappy)"
        con.execute(copy_sql, [dataset_glob, *params])
        print(f"Wrote query output to {output}")
    else:
        print(con.execute(sql, [dataset_glob, *params]).df().head(50).to_string(index=False))


if __name__ == "__main__":
    main()
