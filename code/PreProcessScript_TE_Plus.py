import csv
import os
import sys
import time

import numpy as np
import pandas as pd

from code.path_utils import coerce_input_path, coerce_output_path, locate_file


def load_gene_chr_mapping(filename):
    """Load mapping from gene name to chromosome."""
    path = coerce_input_path(filename)
    with open(path, "r", encoding="utf-8") as f:
        return {
            gene: chromosome
            for chromosome, gene in (line.strip().split() for line in f)
            if not gene.startswith("chr")
        }


def load_file_to_list(filename, delimiter: str = "\n"):
    """Load a file into a list, stripping trailing newlines."""
    path = coerce_input_path(filename)
    with open(path, "r", encoding="utf-8") as f:
        return [line.replace(delimiter, "") for line in f]


def create_TENET_Plus_pairs(Peaks, genes_per_chromosome, gene_names_dict, tf_list):
    """Create gene–gene pairs based on TENET+ logic."""
    gene_pairs = []
    for peak in Peaks:
        peak_chr = peak.split("-")[0]
        peak_idx = gene_names_dict.get(peak, None)
        if peak_idx is not None:
            gene_pairs.extend(
                [
                    (peak_idx + 1, gene_names_dict[gene] + 1)
                    for gene in genes_per_chromosome.get(peak_chr, [])
                    if gene in gene_names_dict
                ]
            )
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend(
            [
                (tf_idx + 1, gene_names_dict[gene] + 1)
                for gene in gene_names_dict
                if gene_names_dict[gene] != tf_idx
            ]
        )
    return gene_pairs


def create_tf_peak_pairs(gene_names_dict, only_peak_names_dict, tf_list):
    """Create TF->peak pairs."""
    peak_pairs = []
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        peak_pairs.extend(
            [
                (tf_idx + 1, only_peak_names_dict[peak] + 1)
                for peak in only_peak_names_dict
                if only_peak_names_dict[peak] != tf_idx
            ]
        )
    return peak_pairs


def create_tf_gene_pairs(only_gene_names_dict, tf_list, gene_names_dict):
    """Create TF->gene pairs."""
    gene_pairs = []
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend(
            [
                (tf_idx + 1, only_gene_names_dict[gene] + 1)
                for gene in only_gene_names_dict
                if only_gene_names_dict[gene] != tf_idx
            ]
        )
    return gene_pairs


def create_tf_gene_peak_pairs(gene_names_dict, tf_list):
    """Create TF->gene and TF->peak pairs in one list."""
    gene_pairs = []
    tf_indices = [gene_names_dict[tf] for tf in tf_list if tf in gene_names_dict]
    for tf_idx in tf_indices:
        gene_pairs.extend(
            [
                (tf_idx + 1, gene_names_dict[gene] + 1)
                for gene in gene_names_dict
                if gene_names_dict[gene] != tf_idx
            ]
        )
    return gene_pairs


def create_cis_peaksource_pairs(Peaks, genes_per_chromosome, gene_names_dict):
    """Create cis peak->gene pairs."""
    gene_pairs = []
    for peak in Peaks:
        peak_chr = peak.split("-")[0]
        peak_idx = gene_names_dict.get(peak, None)
        if peak_idx is not None:
            gene_pairs.extend(
                [
                    (peak_idx + 1, gene_names_dict[gene] + 1)
                    for gene in genes_per_chromosome.get(peak_chr, [])
                    if gene in gene_names_dict
                ]
            )
    return gene_pairs


def load_sif_connections(filename, gene_names_dict):
    """
    Load connections from a SIF file and extract TF->gene, TF->peak,
    and peak->gene pairs.
    """
    tf_gene_pairs = []
    tf_peak_pairs = []
    peak_gene_pairs = []
    path = locate_file(filename)
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline()  # skip header
        for line_number, line in enumerate(f, start=2):
            parts = line.strip().split()
            if len(parts) < 7:
                print(
                    f"Warning: Line {line_number} in SIF file does not have enough "
                    "columns. Skipping."
                )
                continue
            TF, gene, peak, TF_GN_TE, TF_PK_TE, PK_GN_TE, distance = parts
            # TF -> Gene
            if TF in gene_names_dict and gene in gene_names_dict:
                tf_gene_pairs.append(
                    (gene_names_dict[TF] + 1, gene_names_dict[gene] + 1)
                )
            else:
                print(f"Warning: TF or Gene not found at line {line_number}.")
            # TF -> Peak
            if TF in gene_names_dict and peak in gene_names_dict:
                tf_peak_pairs.append(
                    (gene_names_dict[TF] + 1, gene_names_dict[peak] + 1)
                )
            else:
                print(f"Warning: TF or Peak not found at line {line_number}.")
            # Peak -> Gene
            if peak in gene_names_dict and gene in gene_names_dict:
                peak_gene_pairs.append(
                    (gene_names_dict[peak] + 1, gene_names_dict[gene] + 1)
                )
            else:
                print(f"Warning: Peak or Gene not found at line {line_number}.")
    return tf_gene_pairs, tf_peak_pairs, peak_gene_pairs


def create_peak_peak_pairs(Peaks, only_peak_names_dict, max_distance: int | None = 1_000_000):
    """
    Create directed cis peak->peak pairs within an optional distance limit.
    Distance is calculated as the minimal gap between the two peak intervals.
    Both directions (peak1->peak2 and peak2->peak1) are included.
    If max_distance is None, form all cis pairs on the same chromosome.
    """
    peak_pairs = []
    peaks_by_chr = {}
    for peak in Peaks:
        parts = peak.split("-")
        if len(parts) < 3:
            continue
        chrom, start_str, end_str = parts[0], parts[1], parts[2]
        try:
            start, end = int(start_str), int(end_str)
        except ValueError:
            continue
        if peak in only_peak_names_dict:
            peaks_by_chr.setdefault(chrom, []).append((peak, start, end))
    for chrom, plist in peaks_by_chr.items():
        for i in range(len(plist)):
            p1, s1, e1 = plist[i]
            for j in range(len(plist)):
                if i == j:
                    continue
                p2, s2, e2 = plist[j]
                if e1 < s2:
                    dist = s2 - e1
                elif e2 < s1:
                    dist = s1 - e2
                else:
                    dist = 0
                if max_distance is None or dist <= max_distance:
                    idx1 = only_peak_names_dict[p1] + 1
                    idx2 = only_peak_names_dict[p2] + 1
                    peak_pairs.append((idx1, idx2))
    return peak_pairs


def load_and_transpose_matrix(file_name, outputfile_name):
    """Load a Parquet expression matrix, transpose, and save it."""
    source_path = locate_file(file_name)
    output_path = coerce_output_path(outputfile_name)
    expression_data = pd.read_parquet(source_path)
    pd.DataFrame(expression_data.transpose().values).to_parquet(output_path)


def run_preprocess(species: str, mode: str, matrix_file: str, sif_file: str | None) -> None:
    """Run the TENET_Plus preprocessing and write all_pairs.csv."""
    start_time = time.time()

    pair_mode = os.getenv("TENET_PAIR_MODE", "default").strip().lower()

    # For all-pair modes we no longer materialise all_pairs.csv here; the TE core
    # enumerates pairs implicitly in a streaming fashion using gene_names.
    if pair_mode in ("gene_only", "all_feature", "all_pair"):
        print(
            f"[TENET_Plus] pair_mode={pair_mode}: skipping explicit all_pairs.csv "
            "generation; TE core will enumerate full pairs implicitly."
        )
        return

    gene_chr_mapping = load_gene_chr_mapping("gene_chr.txt")

    transposed_matrix_path = coerce_output_path("cell_gene_trsps.parquet")
    if not transposed_matrix_path.exists():
        load_and_transpose_matrix(matrix_file, transposed_matrix_path)

    gene_names = load_file_to_list("gene_names")
    Peaks = load_file_to_list("TE_peak_list.txt")
    tf_list = load_file_to_list(
        f"GO_symbol_{species}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt"
    )

    gene_names_dict = {name: idx for idx, name in enumerate(gene_names)}
    only_peak_names_dict = {
        gene: idx for gene, idx in gene_names_dict.items() if gene in Peaks
    }
    only_gene_names_dict = {
        gene: idx for gene, idx in gene_names_dict.items() if gene not in Peaks
    }

    genes_per_chromosome = {}
    for gene, chromosome in gene_chr_mapping.items():
        if gene in gene_names_dict:
            genes_per_chromosome.setdefault(chromosome, []).append(gene)

    if sif_file:
        print("using Triplet pairs")
        tf_gene_pairs, tf_peak_pairs, peak_gene_pairs = load_sif_connections(
            sif_file, gene_names_dict
        )
        tf_gene_pairs = list(set(tf_gene_pairs))
        tf_peak_pairs = list(set(tf_peak_pairs))
        peak_gene_pairs = list(set(peak_gene_pairs))
        gene_pairs_indices = tf_gene_pairs + tf_peak_pairs + peak_gene_pairs
    else:
        print("Make pairs")
        if mode == "1":
            gene_pairs_indices = create_TENET_Plus_pairs(
                Peaks, genes_per_chromosome, gene_names_dict, tf_list
            )
        elif mode == "2":
            gene_pairs_indices = create_tf_gene_pairs(
                only_gene_names_dict, tf_list, gene_names_dict
            )
        elif mode == "3":
            gene_pairs_indices = create_tf_peak_pairs(
                gene_names_dict, only_peak_names_dict, tf_list
            )
        elif mode == "4":
            gene_pairs_indices = create_tf_gene_peak_pairs(gene_names_dict, tf_list)
        elif mode == "5":
            gene_pairs_indices = create_cis_peaksource_pairs(
                Peaks, genes_per_chromosome, gene_names_dict
            )
        elif mode == "6":
            gene_pairs_indices = create_peak_peak_pairs(
                Peaks, only_peak_names_dict, max_distance=None
            )
        else:
            print("Error: Invalid option for mode.")
            sys.exit(1)

    output_pairs_path = coerce_output_path("all_pairs.csv")
    with open(output_pairs_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(gene_pairs_indices)

    print(f"--- Preprocess time: {time.time() - start_time:.2f} seconds ---")


def main(argv: list[str] | None = None) -> None:
    """Entry point for `python -m code.PreProcessScript_TE_Plus`."""
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) < 3:
        print(
            "Usage: python -m code.PreProcessScript_TE_Plus "
            "<species> <mode> <matrix_file> [sif_file]",
            file=sys.stderr,
        )
        sys.exit(1)
    species, mode, matrix_file, *rest = argv
    sif_file = rest[0] if rest else None
    run_preprocess(species, mode, matrix_file, sif_file)


if __name__ == "__main__":
    main()
