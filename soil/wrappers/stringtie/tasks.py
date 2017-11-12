import pypeliner.commandline as cli


def assemble(bam_file, ref_gtf_file, out_gtf_file, threads=1):
    """ Assemble a transcriptome using stringtie.

    :param bam_file: Path to BAM file with spliced alignments to assemble
    :param ref_gtf_file: Path to GTF file for the reference transcriptome
    :param out_gtf_file: Path where output file will be written in GTF format
    :param threads: Number of threads to use for assembly
    """

    cmd = [
        'stringtie',
        '-o', out_gtf_file,
        '-p', threads,
        '-G', ref_gtf_file,
        bam_file
    ]

    cli.execute(*cmd)
