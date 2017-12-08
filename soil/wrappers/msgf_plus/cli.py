import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-d', '--in-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-s', '--in-mzml-file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out-file', required=True, type=click.Path(resolve_path=True))
@click.option('--split-size', default=int(1e3), type=int)
def search(in_fasta_file, in_mzml_file, out_file, split_size):
    return workflows.create_search_workflow(
        in_fasta_file,
        in_mzml_file,
        out_file,
        split_size=split_size
    )


@click.group(name='msgf-plus')
def msgf_plus():
    pass


msgf_plus.add_command(search)
