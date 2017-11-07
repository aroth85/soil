import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-i', '--bam_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_vcf_file', required=True, type=click.Path(resolve_path=True))
@click.option('-c', '--chromosomes', multiple=True, type=str)
@click.option('-s', '--split_size', default=int(1e7), type=int)
def pileup2snp(bam_file, ref_genome_fasta_file, out_vcf_file, chromosomes, split_size):
    if len(chromosomes) == 0:
        chromosomes = None

    return workflows.create_pileup2snp_workflow(
        bam_file,
        ref_genome_fasta_file,
        out_vcf_file,
        chromosomes=chromosomes,
        split_size=split_size
    )


@click.group()
def varscan():
    pass

varscan.add_command(pileup2snp)
