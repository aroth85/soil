import click
import functools
import pypeliner


def runner(func):
    """ Wrapper function to create a soil runner.
    """
    @functools.wraps(func)
    def func_wrapper(*args, **kwargs):
        config = {
            'maxjobs': kwargs.pop('max_jobs'),
            'nativespec': kwargs.pop('native_spec'),
            'submit': kwargs.pop('submit'),
            'tmpdir': kwargs.pop('working_dir'),
        }

        workflow = func(*args, **kwargs)

        pyp = pypeliner.app.Pypeline(config=config)

        pyp.run(workflow)

    func_wrapper = click.command()(func_wrapper)

    _add_runner_cli_args(func_wrapper)

    return func_wrapper


def _add_runner_cli_args(func):
    """ Add standard pipeline arguments to command line interface for a runner function.

    :param func: Runner function
    """
    click.option('-wd', '--working_dir', required=True, type=click.Path(resolve_path=True))(func)

    click.option('--max_jobs', default=1, type=int)(func)

    click.option('--native_spec', default='', type=str)(func)

    click.option('--submit', default='local', type=click.Choice(['drmaa', 'local']))(func)
