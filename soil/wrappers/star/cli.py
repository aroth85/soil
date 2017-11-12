import click
import os

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-1', '--fastq_file_1', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-2', '--fastq_file_2', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_bam_file', required=True, type=click.Path(resolve_path=True))
@click.option('-t', '--threads', default=1, type=int)
@click.option(
    '-x',
    '--add_xs_tag',
    is_flag=True,
    help='''Set this to support downwstream cufflinks and stringtie analysis.'''
)
def align(fastq_file_1, fastq_file_2, ref_genome_fasta_file, out_bam_file, add_xs_tag, threads):
    ref_genome_fasta_dir = os.path.dirname(ref_genome_fasta_file)

    return workflows.create_align_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_fasta_dir,
        out_bam_file,
        add_xs_tag=add_xs_tag,
        align_threads=threads,
        sort_threads=threads,
    )


@click.group()
def star():
    pass

star.add_command(align)
