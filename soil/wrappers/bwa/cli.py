import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-1', '--fastq_file_1', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-2', '--fastq_file_2', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_bam_file', required=True, type=click.Path(resolve_path=True))
@click.option('-t', '--threads', default=1, type=int)
def mem(fastq_file_1, fastq_file_2, ref_genome_fasta_file, out_bam_file, threads):
    return workflows.create_mem_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_fasta_file,
        out_bam_file,
        align_threads=threads,
        sort_threads=threads
    )


@click.group()
def bwa():
    pass

bwa.add_command(mem)
