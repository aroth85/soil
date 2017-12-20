import os
import pypeliner.commandline as cli
import shutil

import soil.utils.workflow


def get_chrom_variant_file(chrom, in_file, out_file):
    cmd = [
        'bcftools',
        'view',
        '-r', chrom,
        '-O', 'b',
        '-o', out_file
    ]

    cli.execute(*cmd)


def concat_results(in_files, out_file):
    in_files = soil.utils.workflow.flatten_input(in_files)

    tmp_index_files = []

    cmd = ['bcftools', 'concat', '-O', 'z', '-o', out_file]

    for file_name in in_files:
        tmp_index_file_name = file_name + '.csi'

        if not os.path.exists(tmp_index_file_name):
            cli.execute('bcftools', 'index', file_name)

            tmp_index_files.append(tmp_index_file_name)

    cmd.extend(in_files)

    cli.execute(*cmd)

    for file_name in tmp_index_files:
        os.unlink(file_name)


def run_eagle(genetic_map_file, ref_file, target_file, out_file, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_prefix = os.path.join(tmp_dir, 'results')

    if not os.path.exists(ref_file + '.csi'):
        cli.execute('bcftools', 'index', ref_file)

    if not os.path.exists(target_file + '.csi'):
        cli.execute('bcftools', 'index', target_file)

    cmd = [
        'eagle',
        '--geneticMapFile', genetic_map_file,
        '--vcfOutFormat', 'b',
        '--vcfRef', ref_file,
        '--vcfTarget', target_file,
        '--outPrefix', tmp_prefix
    ]

    cli.execute(*cmd)

    shutil.move(tmp_prefix + '.bcf', out_file)

    shutil.rmtree(tmp_dir)
