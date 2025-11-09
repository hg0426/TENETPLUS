from pathlib import Path
import numpy
import sys


def main() -> None:
    input_path = Path(sys.argv[1])
    targets = []
    in_degree = []

    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            temp = line.split()
            if len(temp) < 3:
                continue
            target = temp[2]
            if target not in targets:
                targets.append(target)
                in_degree.append(1)
            else:
                idx = targets.index(target)
                in_degree[idx] += 1

    if not targets:
        print("No edges found in input; nothing to write.")
        return

    sorted_indices = numpy.argsort(in_degree)
    output_path = input_path.with_name(input_path.name + ".indegree.txt")
    with output_path.open("w", encoding="utf-8") as handle:
        for i in range(len(targets)):
            idx = sorted_indices[-i - 1]
            handle.write(f"{targets[idx]}\t{in_degree[idx]}\n")


if __name__ == "__main__":
    main()
