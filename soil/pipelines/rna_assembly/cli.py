import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-b', '--bam_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of RNA-Seq BAM file to use for custom transcriptome prediction. Should have the XS flag.'''
)
@click.option(
    '-g', '--ref_gtf_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of gene annotations in GTF format.'''
)
@click.option(
    '-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path of reference genome in FASTA format.'''
)
@click.option(
    '-o', '--out_prefix', required=True, type=click.Path(resolve_path=True),
    help='''Prefix of output files.'''
)
@click.option(
    '-t', '--threads', default=1, type=int,
    help='''Number of threads to use for parallel steps.'''
)
def rna_assembly(bam_file, ref_genome_fasta_file, ref_gtf_file, out_prefix, threads):
    out_cds_fasta_file = out_prefix + '.cds.fasta'

    out_protein_fasta_file = out_prefix + '.prot.fasta'

    out_gff_file = out_prefix + '.gff3'

    out_gtf_file = out_prefix + '.gtf'

    return workflows.create_assembly_workflow(
        bam_file,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=threads
    )
