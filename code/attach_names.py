
import argparse
import os
import sys
import pandas as pd


def load_names(names_file: str | None, matrix_parquet: str | None, expected_len: int | None = None) -> list[str]:
    names: list[str] | None = None

    # 1) Explicit names file (one name per line)
    if names_file and os.path.exists(names_file):
        with open(names_file, 'r') as f:
            names = [ln.strip() for ln in f if ln.strip()]

    # 2) Try to infer from matrix parquet (index or a dedicated column)
    if names is None and matrix_parquet and os.path.exists(matrix_parquet):
        try:
            dfm = pd.read_parquet(matrix_parquet)
            # Common patterns: gene names as index, or a column named like 'gene' / 'Gene'
            if dfm.index is not None and dfm.index.nlevels == 1 and not isinstance(dfm.index, pd.RangeIndex):
                names = dfm.index.astype(str).tolist()
            else:
                for col in ('gene', 'Gene', 'gene_id', 'GeneID'):
                    if col in dfm.columns:
                        names = dfm[col].astype(str).tolist()
                        break
        except Exception:
            pass

    if names is None:
        raise RuntimeError("Could not load gene names. Provide --names_file or --matrix_parquet with names in index/column.")

    if expected_len is not None and len(names) < expected_len:
        raise RuntimeError(f"Names list too short (len={len(names)}) for expected index up to {expected_len}.")
    return names


def attach_names(df_in: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    # Indices in results are 1-based; convert to 0-based for list access
    try:
        src_idx = (df_in['Source'].astype(int) - 1).to_numpy()
        tgt_idx = (df_in['Target'].astype(int) - 1).to_numpy()
    except KeyError:
        raise RuntimeError("Input file must contain 'Source' and 'Target' columns.")

    max_idx = max(src_idx.max(initial=0), tgt_idx.max(initial=0))
    if max_idx >= len(names):
        raise RuntimeError(f"Names list length {len(names)} is smaller than max index {max_idx+1} in results.")

    df = df_in.copy()
    df['SourceName'] = [names[i] for i in src_idx]
    df['TargetName'] = [names[i] for i in tgt_idx]
    return df


def main():
    ap = argparse.ArgumentParser(description="Attach gene names to TE result files (add SourceName/TargetName columns).")
    ap.add_argument('--in', dest='inp', required=True, help="Input parquet (e.g., TE_result_all.parquet or TE_linear_perm.parquet)")
    ap.add_argument('--out', dest='out', required=False, default=None, help="Output parquet (default: add _with_names suffix)")
    ap.add_argument('--names_file', default='gene_names', help="Text file: one gene name per line in the same order used for indices")
    ap.add_argument('--matrix_parquet', default='cell_gene_trsps.parquet', help="Matrix parquet to infer names (index/column) if names_file not present")
    args = ap.parse_args()

    df = pd.read_parquet(args.inp)
    # Determine the highest index needed (1-based)
    exp_len = int(max(df['Source'].max(), df['Target'].max())) if {'Source','Target'} <= set(df.columns) else None
    names = load_names(args.names_file if os.path.exists(args.names_file) else None, args.matrix_parquet, expected_len=exp_len)
    df2 = attach_names(df, names)

    out = args.out
    if out is None:
        base = args.inp[:-8] if args.inp.endswith('.parquet') else args.inp
        out = f"{base}_with_names.parquet"

    df2.to_parquet(out, index=False)
    print(f"Wrote: {out}  (rows={len(df2):,})")


if __name__ == '__main__':
    main()
