import click
import os

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-f', '--fastq-files', required=True, multiple=True, nargs=2, type=click.Path(exists=True, resolve_path=True),
    help='''A pair of file paths for 1 and 2 file from paired end sequencing. Can be specified multiple times if
    multiple lanes where run.'''
)
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference genome in FASTA format to align against. STAR index files should be in same directory.'''
)
@click.option(
    '-o', '--out-bam-file', required=True, type=click.Path(resolve_path=True),
    help='''Path where output will be written in coordinate sorted duplicate marked BAM format.'''
)
@click.option(
    '-l', '--library-id', default=None, type=str,
    help='''Name of library sequenced to create FASTQ files.'''
)
@click.option(
    '-rg', '--read-group-ids', multiple=True, type=str,
    help='''Read group ID to be used for the lanes. Should be set to match -f, that is the same order and number of
    times. If not set it will be guessed from the file name.'''
)
@click.option(
    '-s', '--sample-id', default=None, type=str,
    help='''Name of sample used to prepare library.'''
)
@click.option(
    '-t', '--threads', default=1, type=int,
    help='''Number of threads used in parallel steps of workflow.'''
)
@click.option(
    '-x', '--add-xs-tag', is_flag=True,
    help='''Add XS tag to BAM alignment. Set this to support downwstream cufflinks and stringtie analysis.'''
)
def align(fastq_files, ref_genome_fasta_file, out_bam_file, add_xs_tag, library_id, read_group_ids, sample_id, threads):
    fastq_files_1, fastq_files_2, read_group_info = soil.utils.cli.parse_alignment_cli_args(
        fastq_files, library_id, read_group_ids, sample_id)

    ref_genome_fasta_dir = os.path.dirname(ref_genome_fasta_file)

    return workflows.create_multiple_lane_align_workflow(
        fastq_files_1,
        fastq_files_2,
        ref_genome_fasta_dir,
        out_bam_file,
        add_xs_tag=add_xs_tag,
        align_threads=threads,
        merge_threads=threads,
        read_group_info=read_group_info,
        sort_threads=threads,
    )


@click.group()
def star():
    pass

star.add_command(align)
