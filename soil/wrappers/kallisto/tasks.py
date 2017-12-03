import os
import pypeliner.commandline as cli
import shutil


def build_index(in_file, out_file, kmer_length=31):
    """ Build an index file for Kallisto
    
        :param in_file: Path to transcriptome fasta file
        :param out_file: Path where Kallisto index file will be written 
    """
    
    cmd = [
        'kallisto',
        'index',
        '-i', out_file,
        '-k', kmer_length,
        in_file,
    ]

    cli.execute(*cmd)


def quantify(
        fastq_file_1,
        fastq_file_2,
        index_file,
        out_file,
        tmp_dir,
        num_bootstraps=100):

    cmd = [
        'kallisto',
        'quant',
        '-i', index_file,
        '-o', tmp_dir,
        '-b', num_bootstraps,
        fastq_file_1,
        fastq_file_2,
    ]

    cli.execute(*cmd)

    tmp_out_file = os.path.join(tmp_dir, 'abundance.tsv')

    shutil.move(tmp_out_file, out_file)

