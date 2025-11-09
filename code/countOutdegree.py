from pathlib import Path
import numpy
import sys


def main() -> None:
    input_path = Path(sys.argv[1])
    tf_list = []
    out_degree = []

    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            temp = line.split()
            if not temp:
                continue
            source = temp[0]
            if source not in tf_list:
                tf_list.append(source)
                out_degree.append(1)
            else:
                idx = tf_list.index(source)
                out_degree[idx] += 1

    if not tf_list:
        print("No edges found in input; nothing to write.")
        return

    sorted_indices = numpy.argsort(out_degree)
    output_path = input_path.with_name(input_path.name + ".outdegree.txt")
    with output_path.open("w", encoding="utf-8") as handle:
        for i in range(len(tf_list)):
            idx = sorted_indices[-i - 1]
            handle.write(f"{tf_list[idx]}\t{out_degree[idx]}\n")


if __name__ == "__main__":
    main()
