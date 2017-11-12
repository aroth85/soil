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
    return workflows.create_assembly_workflow(
        bam_file,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_prefix,
        threads=threads
    )


@click.group()
def rna_assembly():
    pass

rna_assembly.add_command(from_bam)
