import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-n', '--normal_bam_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of BAM file for normal (non-malignant) sample.'''
)
@click.option(
    '-t', '--tumour_bam_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of BAM file for tumour (malignant) sample.'''
)
@click.option(
    '-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of reference genome FASTA file that BAM files were aligned to.'''
)
@click.option(
    '-o', '--out_vcf_file', required=True, type=click.Path(resolve_path=True),
    help='''Path where output file will be written in bgzip compressed VCF format.'''
)
@click.option(
    '-c', '--chromosomes', multiple=True, type=str,
    help='''Chromosome to analyze. Can be specified multiple times i.e. -c chr1 -c chrX to analyze chromosomes 1 and X.
    '''
)
@click.option(
    '-e', '--exome', is_flag=True,
    help='''Set this if the data is from exome sequencing. Disables depth filtering.'''
)
def somatic(normal_bam_file, tumour_bam_file, ref_genome_fasta_file, out_vcf_file, chromosomes, exome):
    if len(chromosomes) == 0:
        chromosomes = None

    return workflows.create_somatic_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_vcf_file,
        chromosomes=chromosomes,
        is_exome=exome,
    )


@click.group()
def strelka():
    pass

strelka.add_command(somatic)
