import click

import soil.utils.cli

import mappability.workflows
import workflows


@soil.utils.cli.runner
@click.option('-o', '--out_dir', required=True, type=click.Path(resolve_path=True))
@click.option('-r', '--ref_genome_version', default='GRCh37', type=click.Choice(['GRCh37', ]))
@click.option('-t', '--threads', default=1, type=int)
@click.option('--cosmic', is_flag=True)
def create(ref_genome_version, out_dir, cosmic, threads):
    """ Download and index reference data.
    """
    return workflows.create_ref_data_workflow(ref_genome_version, out_dir, cosmic=cosmic, threads=threads)


@soil.utils.cli.runner
@click.option('-o', '--out_file', required=True, type=click.Path(resolve_path=True))
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference genome in FASTA format to align against. BWA index files should be in same directory.'''
)
@click.option('-s', '--split-size', default=int(1e7), type=int)
@click.option('-t', '--threads', default=1, type=int)
def mappability(ref_genome_fasta_file, out_file, split_size, threads):
    mappability.workflows.create_mappability_workflow(
        ref_genome_fasta_file,
        out_file,
        k=100,
        max_map_qual=60,
        split_size=split_size,
        threads=threads,
    )
