import pypeliner


def add_pipeline_args(parser):
    """ Add standard pipeline arguments to command line interface.

    :param parser: An argparser.ArgumentParser instance.
    """
    parser.add_argument('-wd', '--working_dir', required=True)

    parser.add_argument('--max_jobs', type=int, default=1)

    parser.add_argument('--native_spec', default='')

    parser.add_argument('--submit', choices=['drmaa', 'local'], default='local')


def run_workflow(args, workflow):
    """ Execute a Pypeliner workflow.

    :param args: Arguments from argparser.ArgumentParser. Should have been run through :py:func:`add_pipeline_args`.
    :param workflow: The Pypeliner workflow to execute.
    """
    config = {
        'maxjobs': args.max_jobs,
        'nativespec': args.native_spec,
        'submit': args.submit,
        'tmpdir': args.working_dir,
    }

    pyp = pypeliner.app.Pypeline(config=config)

    pyp.run(workflow)
