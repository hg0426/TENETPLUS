import argparse
import os
import re
import sys
import numpy as np
import pandas as pd


PEAK_RE = re.compile(r"^chr(?:[0-9]{1,2}|X|Y|M|MT)-([0-9]+)-([0-9]+)$", re.IGNORECASE)


def is_peak(col: str) -> bool:
    return PEAK_RE.match(str(col)) is not None


def load_matrix(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    elif ext == ".csv":
        return pd.read_csv(path, index_col=0)
    elif ext in (".tsv", ".txt"):
        return pd.read_table(path, index_col=0)
    else:
        raise ValueError(f"Unsupported matrix format: {ext}")


def rna_normalise(df: pd.DataFrame, libsize_target: float = 1e4, log1p: bool = True) -> pd.DataFrame:
    # Library size per cell (row sum) over RNA features
    sums = df.sum(axis=1).astype(float)
    # Avoid division by zero
    scale = np.where(sums.to_numpy() > 0, libsize_target / sums.to_numpy(), 0.0)
    df_scaled = df.mul(scale, axis=0)
    if log1p:
        # use np.log1p on values; preserve sparse? DataFrame is dense here
        return np.log1p(df_scaled)
    return df_scaled


def atac_tfidf(df: pd.DataFrame, l2norm: bool = True) -> pd.DataFrame:
    # Term-frequency: divide each cell by its row sum (counts per cell)
    rs = df.sum(axis=1).astype(float)
    tf = df.div(rs.replace(0, np.nan), axis=0).fillna(0.0)
    # Document frequency: number of nonzero cells per feature
    # Use >0 as present
    dfreq = (df > 0).sum(axis=0).astype(float)
    n_docs = df.shape[0]
    # Smooth IDF to avoid div by zero
    idf = np.log1p(n_docs / (1.0 + dfreq))
    tfidf = tf.mul(idf, axis=1)
    if l2norm:
        # L2-normalize per cell
        norms = np.sqrt((tfidf ** 2).sum(axis=1)).replace(0, np.nan)
        tfidf = tfidf.div(norms, axis=0).fillna(0.0)
    return tfidf


def main():
    ap = argparse.ArgumentParser(description="Modality-specific preprocessing for RNA/ATAC and transpose to cell_gene_trsps.parquet")
    ap.add_argument("--input", required=True, help="Input matrix (cells x features) parquet/csv/tsv/txt")
    ap.add_argument("--output", default="cell_gene_trsps.parquet", help="Output transposed matrix path")
    ap.add_argument("--modality", choices=["rna", "atac", "auto"], default="auto", help="Which modality transform to apply. 'auto' applies per-column type")
    ap.add_argument("--rna_libsize_target", type=float, default=1e4, help="Target library size for RNA normalisation")
    ap.add_argument("--no_rna_log1p", action="store_true", help="Disable log1p for RNA")
    ap.add_argument("--no_atac_l2norm", action="store_true", help="Disable L2 normalization for ATAC TFIDF")
    args = ap.parse_args()

    df = load_matrix(args.input)
    # Columns are features; rows are cells
    cols = df.columns

    if args.modality == "rna":
        df_proc = rna_normalise(df, libsize_target=args.rna_libsize_target, log1p=not args.no_rna_log1p)
    elif args.modality == "atac":
        df_proc = atac_tfidf(df, l2norm=not args.no_atac_l2norm)
    else:
        # auto: split by feature type; apply RNA transform to genes; ATAC TFIDF to peaks
        peak_mask = cols.map(lambda c: is_peak(str(c)))
        if peak_mask.any():
            gene_cols = cols[~peak_mask]
            peak_cols = cols[peak_mask]
            df_gene = df[gene_cols]
            df_peak = df[peak_cols]
            df_gene_proc = rna_normalise(df_gene, libsize_target=args.rna_libsize_target, log1p=not args.no_rna_log1p)
            df_peak_proc = atac_tfidf(df_peak, l2norm=not args.no_atac_l2norm)
            df_proc = pd.concat([df_gene_proc, df_peak_proc], axis=1)
            # Preserve original column order
            df_proc = df_proc[cols]
        else:
            # No peaks detected; treat as RNA
            df_proc = rna_normalise(df, libsize_target=args.rna_libsize_target, log1p=not args.no_rna_log1p)

    # Transpose to feature x cell and write parquet
    df_proc.T.to_parquet(args.output)


if __name__ == "__main__":
    main()

