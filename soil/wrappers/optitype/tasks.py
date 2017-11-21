import os
import pypeliner.commandline as cli
import shutil

import soil.utils.file_system


def run_optitype(reads_1_fastq, reads_2_fastq, hla_type_file, tmp_dir, is_rna=False, threads=1):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    config_file = os.path.join(tmp_dir, 'config.ini')

    _write_optitype_config(config_file, threads=threads)

    cmd = [
        'OptiTypePipeline.py',
        '-i', reads_1_fastq, reads_2_fastq,
        '-o', tmp_dir,
        '-v',
        '--config', config_file,
    ]

    if is_rna:
        cmd.append('--rna')

    else:
        cmd.append('--dna')

    cli.execute(*cmd)

    results_subdir = os.listdir(tmp_dir)

    results_subdir = [x for x in results_subdir if not x.endswith('config.ini')]

    assert len(results_subdir) == 1

    date_time_subdir = results_subdir[0]

    results_file = os.path.join(tmp_dir, date_time_subdir, date_time_subdir + '_result.tsv')

    shutil.move(results_file, hla_type_file)

    shutil.rmtree(tmp_dir)


def _write_optitype_config(file_name, threads=1):
    razers_exe = soil.utils.file_system.find('razers3', os.environ['CONDA_PREFIX'])

    from ConfigParser import ConfigParser

    parser = ConfigParser()

    parser.add_section('mapping')
    parser.set('mapping', 'razers3', razers_exe)
    parser.set('mapping', 'threads', threads)

    parser.add_section('ilp')
    parser.set('ilp', 'solver', 'glpk')
    parser.set('ilp', 'threads', 1)

    parser.add_section('behavior')
    parser.set('behavior', 'deletebam', 'true')
    parser.set('behavior', 'unpaired_weight', '0')
    parser.set('behavior', 'use_discordant', 'false')

    with open(file_name, 'w') as fh:
        parser.write(fh)
