import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-1', '--fastq_file_1', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-2', '--fastq_file_2', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_tsv_file', required=True, type=click.Path(resolve_path=True))
@click.option('-t', '--threads', default=1, type=int)
@click.option('--rna', is_flag=True)
def mixcr(fastq_file_1, fastq_file_2, out_tsv_file, threads, rna):
    if rna:
        return workflows.create_rnaseq_workflow(fastq_file_1, fastq_file_2, out_tsv_file, threads=threads)

    else:
        return workflows.create_basic_workflow(fastq_file_1, fastq_file_2, out_tsv_file, threads=threads)
