'''
Created on 23 Mar 2017

@author: Andrew Roth
'''
import os
import pyopenms
import pypeliner.commandline as cli
import re
import shutil


def build_index(in_file, sentinel_file, tda=0):
    cmd = [
        'msgf_plus',
        'edu.ucsd.msjava.msdbsearch.BuildSA',
        '-Xmx4G',
        '-d', in_file,
        '-tda', tda
    ]

    cli.execute(*cmd)

    open(sentinel_file, 'w').close()


def run_search(
        db_file,
        in_file,
        out_file,
        tmp_dir,
        add_decoys=False,
        add_features=False,
        enzyme=1,
        fixed_mods=None,
        fragment_method=0,
        instrument=0,
        isotope_range=(-1, 2),
        max_mods=2,
        num_tolerable_termini=2,
        num_threads=1,
        precursor_mass_tolerance='20ppm',
        variable_mods=None):

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    mod_file = os.path.join(tmp_dir, 'mods.txt')

    tmp_file = os.path.join(tmp_dir, 'search.mzid')

    write_mods_files(mod_file, max_mods, fixed_mods, variable_mods)

    cmd = [
        'msgf_plus',
        '-d', db_file,
        '-s', in_file,
        '-o', tmp_file,
        '-e', enzyme,
        '-inst', instrument,
        '-m', fragment_method,
        '-mod', mod_file,
        '-ntt', num_tolerable_termini,
        '-t', precursor_mass_tolerance,
        '-thread', num_threads,
        '-ti', ','.join([str(x) for x in isotope_range]),
    ]

    if add_decoys:
        cmd.extend(['-tda', 1])

    if add_features:
        cmd.extend(['-addFeatures', 1])

    cli.execute(*cmd)

    shutil.move(tmp_file, out_file)

    shutil.rmtree(tmp_dir)


def write_mods_files(file_name, num_mods, fixed_mods, variable_mods):
    with open(file_name, 'w') as fh:
        fh.write('NumMods={}\n'.format(num_mods))

        if fixed_mods is not None:
            for mod in fixed_mods:
                mod_str = _get_modification_string(mod, fixed=True)

                fh.write('\t'.join((mod_str, '# {}'.format(mod))) + '\n')

        if variable_mods is not None:
            for mod in variable_mods:
                mod_str = _get_modification_string(mod, fixed=False)

                fh.write('\t'.join((mod_str, '# {}'.format(mod))) + '\n')


def _get_modification_string(mod, fixed=True):
    base_mod, aa = re.search('(.*)\s\((.*)\)', mod).groups()

    db = pyopenms.ModificationsDB()

    idx = db.findModificationIndex(mod)

    res_mod = db.getModification(idx)

    weight = res_mod.getDiffFormula().getMonoWeight()

    term_spec = res_mod.getTermSpecificityName(res_mod.getTermSpecificity())

    if term_spec == 'none':
        term_spec = 'any'

    if fixed:
        mod_type = 'fix'

    else:
        mod_type = 'opt'

    return ','.join((
        str(weight),
        aa,
        mod_type,
        term_spec,
        base_mod,
    ))
