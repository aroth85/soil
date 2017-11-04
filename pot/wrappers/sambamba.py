"""
Wrappers for Sambamba http://lomereiter.github.io/sambamba/.
"""
import os
import pypeliner.commandline as cli
import shutil

from pot.utils.workflow import flatten_input


def sort(in_file, out_file, tmp_dir, memory=24, threads=1):
    """ Sort a BAM file by coordinate.

        :param in_file: Path to unsorted BAM file.
        :param out_file: Path where coordinate sorted BAM file will be written.
        :param tmp_dir: Path where a directory will be created to store temporary files.
        :param memory: Rough limit on the amount of memory to use in gigabytes.
        :param threads: Number of threads to use.
    """
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    cmd = [
        'sambamba',
        'sort',
        '-m', '{}G'.format(memory),
        '--tmpdir', tmp_dir,
        '-o', out_file,
        '-t', threads,
        in_file
    ]

    cli.execute(*cmd)

    shutil.rmtree(tmp_dir)


def markdups(in_files, out_file, tmp_dir, threads=1):
    """ Merge files and mark duplicate reads in a file.

        :param in_files: The path of a BAM file or dictionary with values being paths of BAM files. If a dictionary is
            passed the keys are ignored. The BAM files must be coordinate sorted.
        :param out_file: Path where merged and duplicate marked BAM file will be written.
        :param tmp_dir: Path where a directory will be created to store temporary files.
        :param threads: Number of threads to use.
    """
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    cmd = [
        'sambamba',
        'markdup',
        '-t', threads,
        '--tmpdir', tmp_dir,
        '--sort-buffer-size', 20480,
        '--io-buffer-size', 1280,
        '--overflow-list-size', 2000000,
        '--hash-table-size', 2621440
    ]

    cmd.extend(flatten_input(in_files))

    cmd.append(out_file)

    cli.execute(*cmd)

    shutil.rmtree(tmp_dir)
