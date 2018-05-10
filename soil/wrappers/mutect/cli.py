import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-n', '--normal-bam-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of BAM file for normal (non-malignant) sample.'''
)
@click.option(
    '-t', '--tumour-bam-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of BAM file for tumour (malignant) sample.'''
)
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of reference genome FASTA file that BAM files were aligned to.'''
)
@click.option(
    '-o', '--out-vcf-file', required=True, type=click.Path(resolve_path=True),
    help='''Path where output file will be written in bgzip compressed VCF format.'''
)
@click.option(
    '-c', '--chromosomes', multiple=True, type=str,
    help='''Chromosome to analyze. Can be specified multiple times i.e. -c chr1 -c chrX to analyze chromosomes 1 and X.
    '''
)
def paired(normal_bam_file, tumour_bam_file, ref_genome_fasta_file, out_vcf_file, chromosomes, split_size):
    if len(chromosomes) == 0:
        chromosomes = None

    return workflows.create_mutect_paired_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_vcf_file,
        chromosomes=chromosomes,
        normal_name='normal',
        split_size=int(1e7),
        tumour_name='tumour'
    )


@click.group()
def mutect():
    pass


mutect.add_command(mutect)
