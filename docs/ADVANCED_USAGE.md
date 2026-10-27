# Advanced usage

Run these commands from the package root after the main workflow in the
[README](../README.md). Each tool also exposes its arguments with `--help`.

## Configuration files and reruns

Create a config with `./TENET_Plus_for_py.sh --write-config analysis.ini`, or
copy a previous run's `run_config.ini`. The generator refuses to overwrite an
existing file. INI comments start with `#` or `;` on their own line. Values are
literal text; shell substitutions are never executed. Unknown named-section keys and duplicate definitions are rejected.
`[advanced]` accepts supported-prefix environment names directly; spelling and
values there must match the corresponding runtime option. Empty optional entries use runner defaults.

| Section | Settings |
|---|---|
| `[run]` | Required: `matrix`, `jobs`, `trajectory`, `cells`, `history`, `species`, `mode`. Optional: `output` (defaults to `output/` beside the config). |
| `[reference]` | `bundle`, `tf_list`, `feature_types`, `cis_policy` (`legacy`, `human_legacy`, or `strict`). Omit to use the species defaults. |
| `[permutation]` | `mode` (`off`, `grn_edges`, `all_pairs`), `count`, `fdr`, `p_alpha`, `q_alpha`, `candidate_fdr`, `candidate_te_cutoff`. |
| `[localte]` | `enabled`, `chunk_size`, `values_dtype`, `split_export`, `split_output`, `virtual`, `lazy`. Advanced tuning: `buffer_edges`, `export_workers`, `merge_workers`, `read_batch_rows`, `use_threads`. |
| `[network]` | `enabled`, `fdr`, `te_cutoff`, `trim_indirect`, `trim_threshold`, `threads`. |
| `[reuse]` | `mode` (`off`, `on`, `force`), `directory`, `link_mode` (`symlink`, `copy`). |
| `[runtime]` | `python`, `intermediate_save`, `batch_size`, `buffer_rows`. |
| `[advanced]` | Exported runtime/storage environment settings recorded by the runner. Retain these for an unchanged replay; use the named sections for common settings. |

Use `on`/`off` for toggles, a positive integer for `jobs`, and a positive integer
or `auto` for `history`. `mode` is the network mode number 0–9 from the README.
Input paths, the reference bundle/TF/type files, the reuse directory, and the
main output directory resolve relative to the INI file. LocalTE subdirectories
such as `split_output` resolve within the run output directory. An explicit
`--output` is relative to the directory where the command is issued.

For example, edit the generated config to store LocalTE and test GRN candidates:

```ini
[permutation]
mode = grn_edges
count = 100
candidate_fdr = 0.01
candidate_te_cutoff = 0.0
fdr = off
p_alpha = 0.01
q_alpha = 0.05

[localte]
enabled = on
```

`grn_edges` requires an explicit `candidate_fdr` and tests only the TE-derived
GRN candidates at that threshold. `all_pairs` tests all selected pairs. These
are different testing scopes, not interchangeable speed settings. `fdr = on`
selects permutation q-value filtering at `q_alpha`; `off` uses raw p-values at
`p_alpha`. `[network] fdr` separately controls the TE-derived network export.
Neither threshold is silently substituted for the other.

Rerun into a separate directory from any working directory:

```bash
/path/to/TENETPLUS/TENET_Plus_for_py.sh \
  --config /path/to/previous_run/run_config.ini \
  --output /path/to/new_run
```

The snapshot records settings, not copies or content hashes of input data.
Keep input files and references immutable and retain the package/environment
for reproducibility. A rerun starts the workflow again. Prepared-input reuse
below skips preprocessing only; it is not continuation of an interrupted TE
or permutation calculation. `intermediate_save` writes intermediate results
but does not make `--config` a resume command.

### Existing environment-variable settings

The older interactive and seven-argument forms still accept environment
variables. Common mappings are `[run] output` → `TENET_OUTPUT_DIR`,
`[localte] enabled` → `TENET_STORE_LOCAL_TE`, `[network] fdr` → `TENET_GRN_FDR`,
and `[reference] bundle` → `TENET_REFERENCE`. Use the named INI sections for
new batch workflows. `--config` ignores ambient TENET/runtime tuning variables
and applies the file settings, then `--output` if supplied.

## Prepared input reuse

Copy the previous `run_config.ini`, change `[run] output`, and add:

```ini
[reuse]
mode = on
directory = /path/to/previous_run
link_mode = symlink
```

Then run the copied configuration with `--config`. Existing environment-based
scripts may still use `TENET_REUSE_PREPARED`, `TENET_PREPARED_DIR`, and
`TENET_PREPARED_LINK_MODE`.

Reuse checks feature order, species, network mode, reference, TF list, feature
types, and cis policy. Older prepared directories without valid manifests must
be regenerated; `TENET_REUSE_PREPARED=force` does not bypass those checks.
Keep the original inputs immutable: raw matrix, trajectory, and selection
contents are not hashed. If contents change at the same path, use
`TENET_REUSE_PREPARED=off` and a new output directory. Generated preprocessing
files stay in each run's output directory; use distinct directories for
concurrent runs.

## Virtual LocalTE storage

For full edge x timepoint LocalTE analysis, a compact option is the virtual array:

```bash
TENET_OUTPUT_DIR=output/my_run \
TENET_STORE_LOCAL_TE=on \
LOCAL_TE_VIRTUAL=on \
./TENET_Plus_for_py.sh \
  input/merged_expression_data.csv 32 input/trajectory.txt input/cell_select.txt 1 human 1
```

This keeps `TE_result_all.parquet` score-only and writes compact state tables under:

```text
output/my_run/local_te_virtual_array/
output/my_run/local_te_recipe/recipe.json
```

Materialized values preserve the selected value dtype. It skips the reusable
code index by default so the main run stays fast. To materialize conventional
time chunks later, first build the code index and then export with workers:

```bash
python -m code.localte_code_index \
  --virtual_dir output/my_run/local_te_virtual_array \
  --output_dir output/my_run/local_te_code_index

python -m code.localte_virtual_array export \
  --virtual_dir output/my_run/local_te_virtual_array \
  --output_dir output/my_run/local_te_virtual_export \
  --chunk_size 300 \
  --code_index_dir output/my_run/local_te_code_index \
  --workers 16 \
  --merge_workers 8
```

Useful virtual LocalTE variables:

- `LOCAL_TE_VIRTUAL`: `on` or `off` (default `off`)
- `LOCAL_TE_VIRTUAL_OUTPUT_DIR`: output directory (default `local_te_virtual_array`)
- `LOCAL_TE_VIRTUAL_VALUE_DTYPE`: `float32` by default; `float16` reduces precision
- `LOCAL_TE_VIRTUAL_STATE_CODE_BITS`: `32` by default; use `64` only if a run reports state-code overflow
- `LOCAL_TE_VIRTUAL_COMPRESSION`: Parquet compression (default `zstd`)
- `LOCAL_TE_VIRTUAL_BUFFER_EDGES`: virtual row write buffer (default `5000`)
- `LOCAL_TE_VIRTUAL_CODE_INDEX`: `off` by default; turn `on` only when you want the main run to prebuild the export cache


## LocalTE TF activity by group

LocalTE TF activity can be summarized by any timepoint-level grouping, such as
cell type, cluster, condition, lineage, treatment, or developmental stage. The
workflow uses the compact timepoint outdegree table and does not require a
dense edge-by-timepoint LocalTE table.

The optional root `make_GRN_localTE.sh` is an interactive launcher for this
post-processing step; the main TE/permutation workflow does not call it. The
Python command below provides the same underlying network builder for scripts.

First, aggregate TF-GN LocalTE chunks into timepoint outdegree. Direct and
split LocalTE export stores TF-GN chunks in the `TE_TF_GN` subdirectory:

```bash
PYTHONPATH=. \
LOCALTE_OUTPUT_MODE=outdegree \
LOCALTE_TRIM_INDIRECT=yes \
LOCALTE_TRIM_THRESHOLD=0.01 \
LOCAL_TE_OUTDEGREE_HIST=on \
python -m code.make_GRN_localTE \
  --input output/TE_TF_GN.parquet \
  --history_length 1 \
  --fdr 0.01 \
  --te_cutoff 0 \
  --chunk_dir output/local_te_split_chunks/TE_TF_GN \
  --chunk_workers 32 \
  --output_mode outdegree \
  --output_file output/localTE_TF_GN_timepoint_outdegree.parquet
```

Set `LOCALTE_TRIM_INDIRECT=off` to retain untrimmed timepoint outdegree. If a
materialized LocalTE edge table already exists, `trim_indirect_localte_times.py`
can instead trim it and write outdegree with `--outdegree-output` in one pass.

For the activity analysis, provide a CSV, TSV, or parquet map containing
`TimeInt` (or `Time`) and a grouping column. This explicit map is recommended
because it works with arbitrary annotations:

```bash
PYTHONPATH=. \
python -m code.localte_tf_activity_analysis \
  --outdegree output/localTE_TF_GN_timepoint_outdegree.parquet \
  --time-group input/time_to_group.tsv \
  --group-column condition \
  --group-label Condition \
  --outdir output/localte_tf_activity \
  --top-per-group 3 \
  --rolling-window 501 \
  --group-background-alpha 0.24
```

Use a tab-separated mapping covering the timepoints in your outdegree table.
For example:

```text
TimeInt	condition
1	control
2	control
3	treated
```

When group labels are encoded in cell names, they can instead be inferred from
the prepared feature-by-cell parquet with a user-supplied regular expression:

```bash
PYTHONPATH=. \
python -m code.localte_tf_activity_analysis \
  --outdegree output/localTE_TF_GN_timepoint_outdegree.parquet \
  --infer-group-from-cell-gene output/cell_gene_trsps.parquet \
  --history-length 1 \
  --group-regex='^([^_]+)' \
  --group-label Group \
  --outdir output/localte_tf_activity
```

TFs are ranked by activity enrichment inside each group relative to all other
groups. The preferred defaults are `group_enrichment_vs_non_group >= 1.25` and
`group_fraction_of_tf >= 0.15`. If too few TFs pass, fallback selections are
retained and marked with `passed_specificity_filter = False`.

The command writes the time-to-group map, full TF group statistics, selected
group-enriched TFs, smoothed timecourses, a lineplot, and a JSON run summary.


## LocalTE storage and query tools

| Module (`python -m code.<module> --help`) | Use |
|---|---|
| `localte_grn_dataset_query` | Select source, target, and time from a materialized LocalTE network dataset |
| `localte_edge_index_query` | Query the edge-index output and catalog from `make_GRN_localTE --output_mode edge_index` |
| `localte_video_pack` | Pack/unpack LocalTE chunks while preserving their selected dtype |
| `trim_indirect_localte_times` | Trim a materialized timepoint network and optionally write outdegree |

For chunk output, `LOCAL_TE_CHUNK_SIZE=300` and `LOCAL_TE_VALUES_DTYPE=float16`
are the defaults. `LOCAL_TE_MERGE_PARTS=off` keeps worker parts; set it to `on`
only when a single file per chunk is needed. With direct streaming disabled,
export workers default to `min(num_jobs, 16)` and merge workers to
`min(num_jobs, 32)`. `LOCAL_TE_ADVANCED=on` exposes additional interactive tuning.

## Older text TE matrices

The original text-matrix exporter remains available:

```bash
python code/make_GRN_new.py /data/TE_matrix.txt 0.01
```

It reads a headered, tab-separated source-by-target TE matrix, selects positive
scores, and applies a global z-score/normal-tail/BH filter. It is independent
of permutation testing. For current Parquet outputs use Step 4 in the README.

## R post-processing

The optional R helpers are in `code/function_for_tenet.R`. Install their R
packages once, then source the file from the package root:

```r
install.packages(c("ggplot2", "pheatmap", "stringr", "igraph", "dplyr", "jsonlite"))
source("code/function_for_tenet.R")

gene_anno <- read.table("input/gene_anno.bed", header = TRUE)
triplets <- function_PeakSource_Distance(
  Tenet_result_dir = "output",
  gene_anno = gene_anno,
  trim_indirect = TRUE,
  trim_distance = 1000000,
  save = TRUE
)
```

This reads mode-1 networks through `grn_network_outputs.tsv`. `fdr=0.01` must
match the network export threshold. `trim_indirect=TRUE` requires the trimmed
TF->gene network; set it to `FALSE` if trimming was disabled. The older
`TE_result_matrix_row*.fdr0.01.sif` filenames also remain supported.

The human coordinate table is the original GENCODE v44 / GRCh38.p14 file with
62,700 rows. Despite its `.bed` extension, it has a header and these columns:

```text
gene_chr gene_start gene_end gene
chr1 11868 14409 DDX11L2
```

For another species/build, supply the same four-column format with 0-based,
half-open gene intervals and IDs matching the matrix. Do not reuse the human
file. Contig names must match; if aliases differ, pass an explicit named
vector such as `contig_aliases=c("1"="chr1", "MT"="chrM")`. No prefix conversion
is guessed. Where available, peak coordinates and feature types are read from
`feature_manifest.json`; older results use `contig-start-end` or
`contig:start-end` names, including contigs containing hyphens.

Distances are between the **peak interval and gene interval**, not from the
peak to a TSS. Overlapping or adjacent intervals have distance zero. Only
same-contig edges with `distance < trim_distance` are eligible; setting the
cutoff to zero disables the distance cutoff but still excludes different
contigs. Missing gene coordinates and duplicate gene names are excluded with
warnings. Excluding every duplicate gene name preserves the original R
policy, which differs from the TE candidate-generation mapping policy.

The returned columns remain `TF`, `gene`, `peak`, `TF_GN_TE`, `TF_PK_TE`,
`PK_GN_TE`, and `distance`. With `save=TRUE`, the function writes the triplet
`.sif`, a `*_peak_gene_annotations.tsv` table with exclusion reasons, and an
outdegree PDF. Empty input networks yield valid empty outputs. With
`save=FALSE`, results are computed in memory without reading a previous
distance file. The annotation table is also available as
`attr(triplets, "peak_gene_annotations")`.

Older distance results can change for interchromosomal edges: they are now
excluded rather than assigned a distance by comparing unrelated coordinates.
The other original R plotting helpers are retained. The legacy text-matrix
split wrappers `split_mtx_AB` and `split_mtx_C` require an explicit external
`script=` path (or `TENET_LEGACY_SPLIT_AB` / `TENET_LEGACY_SPLIT_C`); current
Parquet splitting is already handled by the main runner.
