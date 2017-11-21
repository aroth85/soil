import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-b', '--bam-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of BAM file to call variants from.'''
)
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of reference genome FASTA file BAM file was aligned to.'''
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
@click.option(
    '--rna', is_flag=True,
    help='''Set this flag if the BAM file contains RNA-Seq data. Will pre-process the files with opossum.'''
)
def single_sample(bam_file, ref_genome_fasta_file, out_vcf_file, chromosomes, rna):
    if len(chromosomes) == 0:
        chromosomes = None

    if rna:
        return workflows.create_rna_single_sample_workflow(
            bam_file,
            ref_genome_fasta_file,
            out_vcf_file,
            chromosomes=chromosomes,
        )

    else:
        return workflows.create_single_sample_workflow(
            bam_file,
            ref_genome_fasta_file,
            out_vcf_file,
            chromosomes=chromosomes,
        )


@click.group()
def platypus():
    pass

platypus.add_command(single_sample)
