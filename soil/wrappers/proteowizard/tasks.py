import os
import pandas as pd
import pypeliner.commandline as cli
import shutil


def split_mzml_file(in_file, out_file_callback, tmp_dir, split_size=int(1e3), validate=False):
    """ Split an mzML file into smaller files of max size `split_size`
    """

    num_specs = _get_num_spectrum(in_file, tmp_dir)

    total = 0

    for i, beg in enumerate(range(0, num_specs, split_size)):
        out_file = out_file_callback[i]

        tmp_file = os.path.splitext(out_file.replace('.tmp', ''))[0]

        tmp_file = tmp_file + '.mzML'

        end = beg + split_size - 1

        cmd = [
            'msconvert',
            in_file,
            '--filter', '"index [{0},{1}]"'.format(beg, end),
            '--outfile', os.path.basename(tmp_file),
            '-o', os.path.dirname(tmp_file),
        ]

        cli.execute(*cmd)

        shutil.move(tmp_file, out_file)

        if validate:
            total += _get_num_spectrum(out_file, tmp_dir)

    if validate:
        assert total == num_specs


def _get_num_spectrum(in_file, tmp_dir):

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    cmd = ['msaccess', in_file, '-x', 'spectrum_table', '-o', tmp_dir]

    cli.execute(*cmd)

    spec_table_file = os.path.basename(in_file) + '.spectrum_table.txt'

    spec_table_file = os.path.join(tmp_dir, spec_table_file)

    spec_table = pd.read_csv(spec_table_file, header=1, sep='\t')

    shutil.rmtree(tmp_dir)

    return spec_table.shape[0]
