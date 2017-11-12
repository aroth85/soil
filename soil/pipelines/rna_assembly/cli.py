import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-b', '--bam_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-g', '--ref_gtf_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_prefix', required=True, type=click.Path(resolve_path=True))
@click.option('-t', '--threads', default=1, type=int)
def from_bam(bam_file, ref_genome_fasta_file, ref_gtf_file, out_prefix, threads):
    out_cds_fasta_file = out_prefix + '.cds.fasta'

    out_protein_fasta_file = out_prefix + '.prot.fasta'

    out_gff_file = out_prefix + '.gff3'

    out_gtf_file = out_prefix + '.gtf'

    return workflows.create_assembly_from_bam_workflow(
        bam_file,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=threads
    )


@soil.utils.cli.runner
@click.option('-1', '--fastq_file_1', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-2', '--fastq_file_2', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-g', '--ref_gtf_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-r', '--ref_genome_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_prefix', required=True, type=click.Path(resolve_path=True))
@click.option('-t', '--threads', default=1, type=int)
def from_fastq(fastq_file_1, fastq_file_2, ref_genome_fasta_file, ref_gtf_file, out_prefix, threads):
    out_bam_file = out_prefix + '.bam'

    out_cds_fasta_file = out_prefix + '.cds.fasta'

    out_protein_fasta_file = out_prefix + '.prot.fasta'

    out_gff_file = out_prefix + '.gff3'

    out_gtf_file = out_prefix + '.gtf'

    return workflows.create_assembly_from_fastq_workflow(
        fastq_file_1,
        fastq_file_1,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_bam_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=threads
    )


@click.group()
def rna_assembly():
    pass

rna_assembly.add_command(from_bam)
rna_assembly.add_command(from_fastq)
