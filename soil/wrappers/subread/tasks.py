import os
import pandas as pd
import pypeliner.commandline as cli
import shutil


def build_gene_counts_file(in_files, out_file):
    """ Merge a list of featureCounts outputs into a gene counts matrix.

        :param in_files: A dicitionary of files to merge keyed on sample name.
        :param out_file: Path where the gene counts matrix file will be written.
    """

    out_df = []

    for sample in in_files:
        df = pd.read_csv(in_files[sample], header=1, sep='\t')

        df = df.drop(['Chr', 'Start', 'End', 'Length', 'Strand'], axis=1)

        df = df.rename(columns={'Geneid': 'gene_id'})

        df = df.set_index('gene_id')

        df.columns = (sample, )

        out_df.append(df)

    out_df = pd.concat(out_df, axis=1)

    out_df.to_csv(out_file, compression='gzip', sep='\t')


def run_feature_counts(in_file, gene_gtf_file, out_file, tmp_dir):
    """ Run featureCounts command to get the number of reads covering genes.

    :param in_file: Path to WTS BAM file to count reads from.
    :param gene_gtf_file: Path to file with gene annotations in GTF format.
    :out_file: Path where output will be written in tsv format.
    :tmp_dir: Temporary director for featureCounts. Will be removed.
    """

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_out_file = os.path.join(tmp_dir, 'counts.tsv')

    cmd = [
        'featureCounts',
        '-a', gene_gtf_file,
        '-o', tmp_out_file,
        in_file
    ]

    cli.execute(*cmd)

    shutil.move(tmp_out_file, out_file)

    shutil.rmtree(tmp_dir)
