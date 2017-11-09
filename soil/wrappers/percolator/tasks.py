'''
Created on 25 Mar 2017

@author: Andrew Roth
'''
import pypeliner.commandline as cli
import os
import shutil


def convert_msgf_to_pin(in_decoy_files, in_target_files, out_file, tmp_space):

    def link_tmp_files(file_list, meta_file, prefix):
        tmp_files = []

        for split_id, file_name in file_list.items():
            tmp_file = os.path.join(tmp_space, '{0}_{1}.mzid'.format(prefix, split_id))

            os.link(file_name, tmp_file)

            tmp_files.append(tmp_file)

        with open(meta_file, 'w') as fh:
            fh.write('\n'.join(tmp_files))

    if os.path.exists(tmp_space):
        shutil.rmtree(tmp_space)

    os.makedirs(tmp_space)

    decoy_meta_file = os.path.join(tmp_space, 'decoy.txt')

    link_tmp_files(in_decoy_files, decoy_meta_file, 'decoy')

    target_meta_file = os.path.join(tmp_space, 'target.txt')

    link_tmp_files(in_target_files, target_meta_file, 'target')

    cmd = [
        'msgf2pin',
        target_meta_file,
        decoy_meta_file,
        '-o', out_file
    ]

    cli.execute(*cmd)

    shutil.rmtree(tmp_space)


def run_percolator(in_file, out_file, db_file=None):
    cmd = [
        'percolator',
        in_file,
        '-X', out_file,
    ]

    if db_file is not None:
        cmd.extend(['-f', db_file])

        cmd.extend(['-P', 'DECOY_'])

    cli.execute(*cmd)
