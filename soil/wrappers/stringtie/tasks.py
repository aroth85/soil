import os
import pypeliner.commandline as cli
import shutil

import soil.utils.workflow


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


def merge(gtf_files, ref_gtf_file, out_gtf_file, tmp_dir, threads=1):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    gtf_list_file = os.path.join(tmp_dir, 'files.txt')

    with open(gtf_list_file, 'w') as fh:
        fh.write('\n'.join(soil.utils.workflow.flatten_input(gtf_files)))

    cmd = [
        'stringtie',
        '--merge',
        '-o', out_gtf_file,
        '-p', threads,
        '-G', ref_gtf_file,
        gtf_list_file
    ]

    cli.execute(*cmd)

    shutil.rmtree(tmp_dir)
