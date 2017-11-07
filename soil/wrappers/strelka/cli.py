import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-n', '--normal_bam_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-t', '--tumour_bam_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_vcf_file', required=True, type=click.Path(resolve_path=True))
@click.option('-c', '--chromosomes', multiple=True, type=str)
@click.option('-s', '--split_size', default=int(1e7), type=int)
def somatic(normal_bam_file, tumour_bam_file, ref_genome_fasta_file, out_vcf_file, chromosomes, split_size):
    if len(chromosomes) == 0:
        chromosomes = None

    return workflows.create_somatic_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_vcf_file,
        chromosomes,
        split_size
    )


@click.group()
def strelka():
    pass

strelka.add_command(somatic)
