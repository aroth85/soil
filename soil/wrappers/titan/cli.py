import click
import os

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
    '-d', '--dbsnp-vcf-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of dbSNP tabix indexed VCF file.'''
)
@click.option(
    '-m', '--mappability-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of bigwig mappability file for reference genome.'''
)
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of reference genome FASTA file that BAM files were aligned to.'''
)
@click.option(
    '-o', '--out-prefix', required=True, type=click.Path(resolve_path=True),
    help='''Path where output files will be written.'''
)
@click.option(
    '-e', '--exome-bed-file', default=None, type=click.Path(resolve_path=True),
    help='''Path of BED file with baits for exome arrays. Must be passed for WES data.'''
)
def titan(
        normal_bam_file,
        tumour_bam_file,
        dbsnp_vcf_file,
        mappability_file,
        ref_genome_fasta_file,
        out_dir,
        exome_bed_file=None):

    out_normal_vcf_file = os.path.join(out_dir, 'normal.vcf.gz')

    out_sentinel_file = os.path.join(out_dir, 'selected', 'done.txt')

    out_stats_file = os.path.join(out_dir, 'stats.tsv')

    return workflows.create_titan_workflow(
        normal_bam_file,
        tumour_bam_file,
        dbsnp_vcf_file,
        mappability_file,
        ref_genome_fasta_file,
        out_sentinel_file,
        out_stats_file,
        exome_bed_file=exome_bed_file,
        out_normal_vcf_file=out_normal_vcf_file
    )
