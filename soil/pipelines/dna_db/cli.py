import click
import os

import soil.ref_data.paths
import soil.utils.cli

import workflows


@soil.utils.cli.runner
@click.option(
    '-n1', '--normal-fastq-file-1', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of normal FASTQ file R1.'''
)
@click.option(
    '-n2', '--normal-fastq-file-2', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of normal FASTQ file R2.'''
)
@click.option(
    '-t1', '--tumour-fastq-file-1', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of tumour FASTQ file R1.'''
)
@click.option(
    '-t2', '--tumour-fastq-file-2', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of tumour FASTQ file R2.'''
)
@click.option(
    '-r', '--ref-data-dir', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of SOIL reference data dir.'''
)
@click.option(
    '-o', '--out-dir', required=True, type=click.Path(resolve_path=True),
    help='''Path of directory where output files will be written.'''
)
@click.option(
    '-e', '--is-exome', is_flag=True,
    help='''Set this if the data is from exome sequencing. Disables depth filtering for Strelka.'''
)
@click.option(
    '-t', '--threads', default=1, type=int,
    help='''Number of threads to use for parallel steps.'''
)
def from_fastq(
        normal_fastq_file_1,
        normal_fastq_file_2,
        tumour_fastq_file_1,
        tumour_fastq_file_2,
        ref_data_dir,
        out_dir,
        is_exome=False,
        threads=1):

    custom_proteome_file = os.path.join(out_dir, 'db.fasta')

    normal_bam_file = os.path.join(out_dir, 'normal.bam')

    tumour_bam_file = os.path.join(out_dir, 'tumour.bam')

    platypus_file = os.path.join(out_dir, 'platypus.vcf.gz')

    strelka_file = os.path.join(out_dir, 'strelka.vcf.gz')

    ref_data_paths = soil.ref_data.paths.SoilRefDataPaths(ref_data_dir)

    return workflows.create_custom_dna_proteome_from_fastq_workflow(
        normal_fastq_file_1,
        normal_fastq_file_2,
        tumour_fastq_file_1,
        tumour_fastq_file_2,
        ref_data_paths.bwa_genome_fasta_file,
        ref_data_paths.proteome_fasta_file,
        normal_bam_file,
        tumour_bam_file,
        custom_proteome_file,
        platypus_file,
        strelka_file,
        genome_version=ref_data_paths.config['pyensembl']['genome'],
        is_exome=is_exome,
        pyensembl_cache_dir=ref_data_paths.pyensembl_cache_dir,
        threads=threads
    )


@click.group(name='dna-db')
def dna_db():
    pass


dna_db.add_command(from_fastq)
