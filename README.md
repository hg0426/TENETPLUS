# TENETPLUS

TENETPLUS reconstructs Transfer Entropy (TE)-based causal gene networks from pseudotime-ordered single-cell transcriptomic and epigenetic data, with an emphasis on RNA + ATAC multiome assays (the non-multiome TENET_TF mode is still supported).

## Method

![tenetplus_workflow](https://github.com/user-attachments/assets/36e6c47a-4f12-4a6a-a74e-1da6eea6ab32)

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014  
https://github.com/neocaleb/TENET

## Highlights
- **Multiome-aware pipeline**: integrates TF->gene, TF->peak, and peak->gene signals while keeping the original TENET_TF flow available.
- **Modern interface**: `TENET_Plus_for_py.sh` uses interactive prompts but also accepts a full CLI to keep batch workflows reproducible.
- **Output consistency**: downstream scripts in `code/` build networks (`make_GRN_new.py`), trim indirect edges, and summarize node degrees.

## Requirements
- `python3`
- `openmpi` (>= 4.0 for MPI-backed TE screening)
- `numpy`, `pandas`, `scipy` (available via PyPI)
- POSIX shell (the wrapper script lives at the repo root)

## Directory layout
- `code/`: Python utilities and the main TE estimation helpers.
- `input/`: place CSV/TSV/RParquet expression matrices, trajectory lists, and cell-select masks here before running TENETPLUS.
- `output/`: default target for TE matrices and GRN files. The runner creates the directory if it does not already exist.

## Step 1 - run `TENET_Plus_for_py.sh`

The shell wrapper is the entry point for both multiome and TF-only modes. Launch it without arguments to walk through the guided prompt sequence, or supply the arguments below to run non-interactively.

### Interactive mode

```bash
./TENET_Plus_for_py.sh
```

Prompts ask for the input matrix, number of parallel jobs, trajectory file, cell-select mask, history length (`k`), species, mode code, and advanced knobs such as modality preprocessing, screening estimator, refinement, permutation testing, local TE export, and time subsampling.

### Batch (CLI) usage

```bash
./TENET_Plus_for_py.sh <input_matrix> <num_jobs> <trajectory_file> <cell_select_file> <history_k> <species> <mode_code> [modality] [screen_mode] [refine_method] [refine_topk] [refine_top_pct] [permute] [perm_n] [perm_topk] [perm_top_pct] [perm_fdr] [perm_q_alpha] [perm_alpha] [store_local_te] [time_stride] [time_pct] [time_seed]
```

- `history_k`: number of previous timepoints to include (commonly 1 or 2).
- `mode_code`: see the table below.
- Optional flags (default in parentheses): `modality` (`auto` for Plus modes, `rna` for TENET_TF; `none` skips preprocessing), `screen_mode` (`kernel`), `refine_method` (`none`), `refine_topk` (0), `refine_top_pct` (0), `permute` (`off`), `perm_n` (100), `perm_topk` (0 = all), `perm_top_pct` (0), `perm_fdr` (`off`), `perm_q_alpha` (0.05), `perm_alpha` (0.01), `store_local_te` (`off`), `time_stride` (1), `time_pct` (100), `time_seed` (42).

### Mode codes

```
0 - TENET_TF: RNA-only pipeline.
1 - TENET_Plus: RNA + ATAC, full TE matrix (TF->gene + TF->peak + peak source).
2 - TENET_Plus rowTF_colGN: only TF -> gene edges.
3 - TENET_Plus rowTF_colPK: only TF -> peak edges.
4 - TENET_Plus rowTF_colGN+PK: both TF -> gene and TF -> peak.
5 - TENET_Plus rowPeak (cis peak-source): builds peak->gene cis interactions.
6 - TENET_Plus peak->peak (cis): builds local peak->peak networks.
```

### Environment tunables

Set these environment variables (export before running the script) to tweak local-TE checkpointing and chunking:

- `LOCAL_TE_CHUNK_SIZE` (default 300): number of timepoints processed together.
- `LOCAL_TE_BUFFER_EDGES` (default 2000): padding per chunk to capture edge effects.
- `LOCAL_TE_EXPORT_WORKERS` (default 0/`LOCAL_TE_MERGE_WORKERS` is 0): controls multiprocessing during export/merge.
- `LOCAL_TE_READ_BATCH_ROWS` (default 8192): rows per read call when streaming.
- `LOCAL_TE_USE_THREADS` (`on`/`off`, default `on`): switch per-chunk threading.
- `LOCAL_TE_VALUES_DTYPE` (default `float16`): dtype for local TE matrices.
- `LOCAL_TE_MERGE_WORKERS` (default empty ≈ 0): worker count for merging local TE files.
- `LOCAL_TE_SPLIT_EXPORT` (`on`): keep per-selector chunk directories (useful when supplying cell-select files).
- `LOCAL_TE_SPLIT_OUTPUT_DIR` (default `local_te_split_chunks`): destination for split exports.

### Output files

All TE matrix outputs are placed under `output/` (or the path you provide via `OUTPUT_DIR`). The full run produces:

- `TE_result_matrix_rowTF_colPK.txt`
- `TE_result_matrix_rowTF_colGN.txt`
- `TE_result_matrix_rowPeak_colGN.txt`
- `TE_result_matrix.txt` (summary of all modes)

Each file can feed into the GRN reconstruction steps below.

## Step 2 - reconstruct the GRN

Use `code/make_GRN_new.py` to convert TE matrices into graph edge lists:

```bash
python code/make_GRN_new.py TE_result_matrix_rowTF_colGN.txt 0.01
python code/make_GRN_new.py TE_result_matrix_rowTF_colPK.txt 0.01
```

The second argument is the TE cutoff (e.g., 0.01). Outputs have names like `TE_result_matrix_rowTF_colGN.fdr0.01.sif`.

## Step 3 - trim indirect edges (TENET original step)

```bash
python code/trim_indirect.py TE_result_matrix_rowTF_colGN.fdr0.01.sif -0.01
```

`trim_indirect.py` removes redundant edges. The cutoff is typically between -0.1 and 0.1; use the same value you found suitable in previous TENET releases.

## Step 4 - summarize degrees

```bash
python code/countOutdegree.py TE_result_matrix_rowTF_colPK.fdr0.01.sif
python code/countOutdegree.py TE_result_matrix_rowTF_colGN.fdr0.01.trimIndirect-0.01.sif
```

Each run writes `<network>.outdegree.txt`.

