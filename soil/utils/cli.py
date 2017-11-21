import click
import functools
import os
import pypeliner
import shutil


def runner(func):
    """ Wrapper function to create a soil runner.
    """
    @functools.wraps(func)
    def func_wrapper(*args, **kwargs):
        resume = kwargs.pop('resume')

        working_dir = kwargs.pop('working_dir')

        if os.path.exists(working_dir) and (not resume):
            raise Exception('''
                Runner failing because working directory {} exists.
                Either remove working directory or use --resume flag to resume an interrupted run.
                '''.format(working_dir))

        config = {
            'maxjobs': kwargs.pop('max_jobs'),
            'nativespec': kwargs.pop('native_spec'),
            'submit': kwargs.pop('submit'),
            'tmpdir': working_dir,
        }

        workflow = func(*args, **kwargs)

        pyp = pypeliner.app.Pypeline(config=config)

        pyp.run(workflow)

        shutil.rmtree(working_dir)

    func_wrapper = click.command(context_settings={'max_content_width': 120})(func_wrapper)

    _add_runner_cli_args(func_wrapper)

    return func_wrapper


def _add_runner_cli_args(func):
    """ Add standard pipeline arguments to command line interface for a runner function.

    :param func: Runner function
    """
    click.option(
        '-wd', '--working_dir', required=True, type=click.Path(resolve_path=True),
        help='''Working directory for runner. Analysis will fail if it exists unless --resume flag is used. Will be
         deleted when analysis is finished.
        '''
    )(func)

    click.option(
        '--max_jobs', default=1, type=int,
        help='''Maximum number of jobs to run.'''
    )(func)

    click.option(
        '--native_spec', default='', type=str,
        help='''String specifying cluster submission parameters. Special values are {mem} for memory requests and
         {threads} for thread requests.

        Examples:
            For single threaded workflows

            "-cwd -V -q byslot.q -l mem_free={mem}G,h_vmem={mem}G"

            For multi-threaded workflows

            "-cwd -V -pe smp {threads} -l mem_free={mem}G,h_vmem={mem}G"
        '''
    )(func)

    click.option(
        '--submit', default='local', type=click.Choice(['drmaa', 'local']),
        help='''Job submission strategy. Use local to run on host machine or drmaa to submit to a cluster.'''
    )(func)

    click.option(
        '--resume', is_flag=True,
        help='''Set this flag if an analysis was interrupted and you would like to resume. Only has an effect if the
         working directory exists.
        '''
    )(func)


def parse_alignment_cli_args(fastq_files, library_id=None, read_group_ids=None, sample_id=None):
    """ Parse standard alignment arguments. Assumes paired end data which all comes from the same library.

    :param fastq_files: List of tuples where each tuple is one FASTA file of a paired end library.
    :param library_id: Identifier of library sequenced to generate FASTA files.
    :param read_group_ids: List of read group identifiers. Length should match length of fastq_files.
    :param sample_id: Identifier of sample used to prepare library.
    """
    fastq_files_1 = {}

    fastq_files_2 = {}

    read_group_info = {}

    if library_id is None:
        library_id = 'UKNOWN'

    if sample_id is None:
        sample_id = 'UKNOWN'

    for idx, (f1, f2) in enumerate(fastq_files):
        if len(read_group_ids) > 0:
            key = read_group_ids[idx]

        else:
            key = os.path.basename(f1).split('.')[0]

        fastq_files_1[key] = f1

        fastq_files_2[key] = f2

        read_group_info[key] = {'ID': key, 'LB': library_id, 'SM': sample_id}

    return fastq_files_1, fastq_files_2, read_group_info
