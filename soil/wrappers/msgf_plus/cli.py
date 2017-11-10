import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-d', '--in_fasta_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-s', '--in_mzml_file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out_file', required=True, type=click.Path(resolve_path=True))
@click.option('--split_size', default=int(1e3), type=int)
def search(in_fasta_file, in_mzml_file, out_file, split_size):
    return workflows.create_search_workflow(
        in_fasta_file,
        in_mzml_file,
        out_file,
        split_size=split_size
    )


@click.group()
def msgf_plus():
    pass

msgf_plus.add_command(search)
