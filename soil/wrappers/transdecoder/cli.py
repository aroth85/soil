import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-g', '--gtf-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to GTF file created by Stringtie of Cufflinks.'''
)
@click.option(
    '-G', '--ref-gtf-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference GTF file.'''
)
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference genome FASTA file.'''
)
@click.option(
    '-o', '--out-prefix', required=True, type=click.Path(resolve_path=True),
    help='''Prefix of output files to generate.'''
)
def transdecoder(gtf_file, ref_gtf_file, ref_genome_fasta_file, out_prefix):
    out_alignment_gff_file = out_prefix + '.gff'

    out_cdna_fasta_file = out_prefix + '.cdna.fasta'

    out_cds_fasta_file = out_prefix + '.cds.fasta'

    out_protein_fasta_file = out_prefix + '.prot.fasta'

    return workflows.create_transdecoder_workflow(
        gtf_file,
        ref_gtf_file,
        ref_genome_fasta_file,
        out_alignment_gff_file,
        out_cdna_fasta_file,
        out_cds_fasta_file,
        out_protein_fasta_file
    )
