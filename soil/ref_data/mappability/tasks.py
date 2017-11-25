from __future__ import division

from Bio import SeqIO
from collections import defaultdict

import pandas as pd
import pypeliner.commandline as cli
import pysam

import soil.utils.workflow


def split_fasta_by_chrom(in_file, out_file_callback):

    with open(in_file, 'r') as in_fh:
        for record in SeqIO.parse(in_fh, format='fasta'):
            out_file = out_file_callback[record.id]

            with open(out_file, 'w') as out_fh:
                SeqIO.write(record, out_fh, format='fasta')


def create_kmer_reads(in_file, out_file_callback, k=100, split_size=int(1e6)):

    with open(in_file, 'r') as in_fh:
        record = SeqIO.parse(in_fh, format='fasta').next()

    chrom = record.id

    seq = str(record.seq)

    file_idx = 0

    file_size = 0

    out_fh = open(out_file_callback[file_idx], 'w')

    for beg in range(len(record.seq) - k):
        if file_size >= split_size:
            out_fh.close()

            file_idx += 1

            file_size = 0

            out_fh = open(out_file_callback[file_idx], 'w')

        end = beg + k

        kmer = seq[beg:end].upper()

        assert len(kmer) == k

        out_fh.write('>{chrom}:{beg}-{end}\n'.format(**locals()))

        out_fh.write(kmer + '\n')

        file_size += 1

    out_fh.close()


def bwa_mem_align(in_file, ref_genome_fasta_file, out_file, threads=1):

    cmd = [
        'bwa', 'mem',
        '-t', threads,
        ref_genome_fasta_file, in_file,
        '|',
        'samtools', 'view', '-b', '-',
        '>',
        out_file
    ]

    cli.execute(*cmd)


def compute_mappability(in_file, out_file, max_map_qual=None):

    bam = pysam.AlignmentFile(in_file)

    probs = defaultdict(float)

    count = defaultdict(int)

    if max_map_qual is None:
        norm_const = 1
    else:
        norm_const = (1 - 10 ** (-max_map_qual / 10))

    for read in bam:
        chrom, region = read.query_name.split(':')

        beg, end = region.split('-')

        beg = int(beg)

        end = int(end)

        q = read.mapping_quality

        p = 1 - 10 ** (-q / 10)

        p = p / norm_const

        for coord in range(beg, end):
            probs[(chrom, coord)] += p

            count[(chrom, coord)] += 1

    probs = pd.Series(probs)

    count = pd.Series(count)

    df = pd.concat([probs, count], axis=1)

    df = df.reset_index()

    df.columns = 'chrom', 'coord', 'mappability', 'count'

    df.to_csv(out_file, index=False, sep='\t')


def compute_chrom_mean_mappability(in_files, norm_const, out_file):
    def collapse_seg(df):
        return pd.Series(
            data=[df['chrom'].iloc[0], df['coord'].min(), df['coord'].max() + 1, df['mappability'].iloc[0]],
            index=['chrom', 'beg', 'end', 'mappability'],
        )

    data = []

    for file_name in soil.utils.workflow.flatten_input(in_files):
        data.append(pd.read_csv(file_name, sep='\t'))

    data = pd.concat(data)

    data = data.groupby(['chrom', 'coord']).sum()

    data['mappability'] = data['mappability'] / data['count']

    data = data.drop('count', axis=1)

    data = data.round(decimals=2)

    data = data.reset_index()

    data.sort_values(by='coord', inplace=True)

    data['seg'] = (data['mappability'].diff() != 0).cumsum()

    data = data.groupby('seg').apply(collapse_seg)

    data.to_csv(out_file, index=False, sep='\t')


def write_bed(in_files, out_file):
    for file_name in soil.utils.workflow.flatten_input(in_files):
        df = pd.read_csv(file_name, sep='\t')

        df.to_csv(out_file, header=False, index=False, mode='a', sep='\t')


def write_chrom_sizes(in_file, out_file):
    sizes = []

    with open(in_file, 'r') as in_fh:
        for record in SeqIO.parse(in_fh, format='fasta'):
            sizes.append({'chrom': record.id, 'size': len(record.seq)})

    sizes = pd.DataFrame(sizes, columns=['chrom', 'size'])

    sizes.to_csv(out_file, header=False, index=False, sep='\t')
