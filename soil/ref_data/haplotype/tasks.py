import pandas as pd


def write_chrom_map_file(in_file, out_file):
    df = pd.read_csv(in_file, sep='\t')

    df = df[['ncbi', 'gencode']]

    df.to_csv(out_file, index=False, header=False, sep='\t')
