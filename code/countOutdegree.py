#!/usr/bin/env python3
"""Count outgoing edges in a TENET SIF network file.

The expected SIF layout is the TENET convention:

    Source<TAB>TE<TAB>Target

The script writes `<input>.outdegree.txt` with two columns: node and outdegree.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Count source outdegree from a TENET SIF file.")
    parser.add_argument("sif", type=Path, help="Input SIF file: Source<TAB>TE<TAB>Target.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path. Defaults to <sif>.outdegree.txt.",
    )
    args = parser.parse_args()

    counts: Counter[str] = Counter()
    with args.sif.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            source = parts[0]
            if source:
                counts[source] += 1

    output = args.output or args.sif.with_suffix(args.sif.suffix + ".outdegree.txt")
    with output.open("w") as out:
        out.write("Node\tOutdegree\n")
        for node, degree in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            out.write(f"{node}\t{degree}\n")

    print(f"Wrote outdegree table to {output}")


if __name__ == "__main__":
    main()
