"""
Wrappers for BWA http://bio-bwa.sourceforge.net/.
"""
import pypeliner.commandline as cli


def index(ref_genome_fasta_file, out_sentinel_file):
    """ Build an index of a FASTA file for use with the BWA alignment programs.

        :param ref_genome_fasta_file: Path to FASTA file for a reference genome which will be aligned against.
        :param out_sentinel_file: Path of sentinel file to show indexing was successful.
    """
    cmd = ['bwa', 'index', ref_genome_fasta_file]

    cli.execute(*cmd)

    open(out_sentinel_file, 'w').close()


def mem_paired_end(fastq_file_1, fastq_file_2, ref_genome_fasta_file, out_bam_file, read_group_info=None, threads=1):
    """ Align paired end FASTQ files using `bwa mem`.

        :param fastq_file_1: Path to FASTQ (or FASTQ.gz) file with 1st paired end reads.
        :param fastq_file_2: Path to FASTQ (or FASTQ.gz) file with 2nd paired end reads.
        :param ref_genome_fasta_file: Path to FASTA file for the reference to be aligned against. Must have BWA index
            files in the same directory.
        :param out_bam_file: Path where output file will be written in BAM format.
        :param readgroup_info: Dictionary containing information about the readgroup. Must have an ID entry.
        :param threads: Number of threads BWA will use.
    """

    cmd = [
        'bwa', 'mem',
        '-t', threads,
    ]

    if read_group_info is not None:
        read_group_str = ['@RG', 'ID:{0}'.format(read_group_info['ID'])]

        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue

            read_group_str.append(':'.join((key, value)))

        read_group_str = '\t'.join(read_group_str)

        cmd.extend(['-R', read_group_str])

    cmd.extend([
        ref_genome_fasta_file, fastq_file_1, fastq_file_2,
        '|',
        'samtools', 'view', '-b', '-',
        '>',
        out_bam_file
    ])

    cli.execute(*cmd)
