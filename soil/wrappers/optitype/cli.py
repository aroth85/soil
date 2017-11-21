import click

import soil.utils.cli
import workflows


@soil.utils.cli.runner
@click.option(
    '-b', '--bam_file', required=True, type=click.Path(exists=True, resolve_path=True),
    help='''Path to BAM file to analyze. Should be a normal (non-malignant) sample if available.'''
)
@click.option(
    '-o', '--out_file', required=True, type=click.Path(resolve_path=True),
    help='''Path where output with HLA information will be written.'''
)
@click.option(
    '--rna', is_flag=True,
    help='''Set this flag if the data is from RNA (WTS). Otherwise assumed to be DNA (WGS/WES).'''
)
@click.option(
    '-th', '--threads', default=1, type=int,
    help='''Number of threads used in parallel steps of workflow.'''
)
def optitype(bam_file, out_file, rna, threads):
    return workflows.create_optitype_workflow(bam_file, out_file, is_rna=rna, threads=threads)
