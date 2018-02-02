import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option('-d', '--in-fasta-file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-s', '--in-mzml-file', required=True, type=click.Path(exists=True, resolve_path=True))
@click.option('-o', '--out-file', required=True, type=click.Path(resolve_path=True))
@click.option('-fm', '--fixed-mods', multiple=True)
@click.option('-vm', '--variable-mods', multiple=True)
@click.option('--max-mods', default=1, type=int)
@click.option('--percolator', is_flag=True)
@click.option('--split-size', default=int(1e3), type=int)
def search(in_fasta_file, in_mzml_file, out_file, fixed_mods, max_mods, percolator, split_size, variable_mods):
    if percolator:
        return workflows.create_percolator_workflow(
            in_fasta_file,
            in_mzml_file,
            out_file,
            fixed_mods=fixed_mods,
            max_mods=max_mods,
            split_size=split_size,
            variable_mods=variable_mods
        )

    else:
        return workflows.create_search_workflow(
            in_fasta_file,
            in_mzml_file,
            out_file,
            fixed_mods=fixed_mods,
            max_mods=max_mods,
            split_size=split_size,
            variable_mods=variable_mods
        )


@click.group(name='msgf-plus')
def msgf_plus():
    pass


msgf_plus.add_command(search)
