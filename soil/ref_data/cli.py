import click
import pkg_resources
import yaml

import soil.ref_data.mappability.workflows
import soil.ref_data.workflows
import soil.utils.cli


@soil.utils.cli.runner
@click.option('-o', '--out_dir', required=True, type=click.Path(resolve_path=True))
@click.option('-c', '--config_file', default=None)
@click.option('-r', '--ref_genome_version', default='GRCh37', type=click.Choice(['GRCh37', ]))
@click.option('-t', '--threads', default=1, type=int)
@click.option('--cosmic', is_flag=True)
def create(config_file, ref_genome_version, out_dir, cosmic, threads):
    """ Download and index reference data.
    """
    if config_file is None:
        config_file = _load_config_file(ref_genome_version)

    with open(config_file, 'r') as fh:
        config = yaml.load(fh)

    return soil.ref_data.workflows.create_ref_data_workflow(config, out_dir, cosmic=cosmic, threads=threads)


@soil.utils.cli.runner
@click.option('-o', '--out_file', required=True, type=click.Path(resolve_path=True))
@click.option(
    '-r', '--ref-genome-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference genome in FASTA format to align against. BWA index files should be in same directory.'''
)
@click.option('-s', '--split-size', default=int(1e7), type=int)
@click.option('-t', '--threads', default=1, type=int)
def mappability(ref_genome_fasta_file, out_file, split_size, threads):
    return soil.ref_data.mappability.workflows.create_mappability_workflow(
        ref_genome_fasta_file,
        out_file,
        k=100,
        max_map_qual=60,
        split_size=split_size,
        threads=threads,
    )


@click.command()
@click.option('-r', '--ref-genome-version', default='GRCh37', type=click.Choice(['GRCh37', ]))
def show_config(ref_genome_version):
    """ Show the YAML config for a reference version.
    """
    config_file = _load_config_file(ref_genome_version)

    with open(config_file, 'r') as fh:
        print(fh.read())


def _load_config_file(ref_genome_version):
    return pkg_resources.resource_filename(
        'soil',
        'ref_data/configs/{}.yaml'.format(ref_genome_version)
    )
