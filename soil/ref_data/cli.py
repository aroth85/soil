import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-o', '--out_dir', required=True, type=click.Path(resolve_path=True))
@click.option('-r', '--ref_genome_version', default='GRCh37', type=click.Choice(['GRCh37', ]))
@click.option('--cosmic', is_flag=True)
def create(ref_genome_version, out_dir, cosmic):
    """ Download and index reference data.
    """
    return workflows.create_ref_data_workflow(ref_genome_version, out_dir, cosmic=cosmic)
