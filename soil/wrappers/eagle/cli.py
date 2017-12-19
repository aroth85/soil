import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-g', '--genetic-map-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path EAGLE genetic map file.'''
)
@click.option(
    '-r', '--ref-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to reference panel in BCF or VCF format.'''
)
@click.option(
    '-i', '--target-file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to target file in BCF or VCF format.'''
)
@click.option(
    '-o', '--out-file', required=True, type=click.Path(resolve_path=True),
    help='''Path where output will be written in BCF format with an index.'''
)
def eagle(genetic_map_file, ref_file, target_file, out_file):
    return workflows.create_ref_panel_phase_workflow(genetic_map_file, ref_file, target_file, out_file)
