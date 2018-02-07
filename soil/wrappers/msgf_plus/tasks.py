import os
import pandas as pd
import pyteomics.fasta
import pypeliner.commandline as cli
import re
import shutil

import soil.utils.package_data
import soil.utils.workflow


def build_decoy_db(in_file, out_file, decoy_only=True, decoy_prefix='XXX_'):
    pyteomics.fasta.write_decoy_db(in_file, out_file, decoy_only=decoy_only, prefix=decoy_prefix)


def _build_index(in_file, tda=0):
    os.environ['MALLOC_ARENA_MAX'] = '2'

    cmd = [
        'msgf_plus',
        'edu.ucsd.msjava.msdbsearch.BuildSA',
        '-Xmx4G',
        '-d', in_file,
        '-tda', tda
    ]

    cli.execute(*cmd)


def build_index_sentinel(in_file, sentinel_file, **kwargs):
    out_file = sentinel_file.replace('.sentinel.tmp', '.fasta')

    build_index(in_file, out_file, **kwargs)

    open(sentinel_file, 'w').close()


def build_index(in_file, out_file, add_decoys=True):
    """ Build an indexed database.

    When using MSGF+ out_file can be set as the search database.

    Parameters
    ----------
    in_file: str
        Path to FASTA format database file.
    out_file: str
        Path where database will be written in FASTA with MSGF+ index files alongside.
    add_decoys: bool
        Whether to add decoy sequences to the database.
    """
    tmp_file = out_file.replace('.tmp', '')

    if add_decoys:
        build_decoy_db(in_file, tmp_file, decoy_only=False, decoy_prefix='XXX_')

    else:
        shutil.copy(in_file, tmp_file)

    _build_index(tmp_file, tda=0)

    shutil.move(tmp_file, out_file)


def clean_up(db_sentinel_files, in_file, out_file):
    """ Clean up the index files for MSGF+.

    Parameters
    ----------
    db_file: str,iterable
        Path(s) to db files used by MSGF+ workflow.
    in_file: str
        Path to final output file of workflow. Should be a tmp.
    out_file: str
        Path to final output of worklow. Should be the real target.
    """
    shutil.copy(in_file, out_file)

    index_exts = ['csarr', 'canno', 'cseq', 'cnlcp']

    for sentinel_file in soil.utils.workflow.flatten_input(db_sentinel_files):
        file_name = sentinel_file.replace('.sentinel', '.fasta')

        for ext in index_exts:
            index_file = file_name.replace('.fasta', '.{}'.format(ext))

            os.unlink(index_file)


def convert_mzid_to_tsv(in_file, out_file):
    tmp_file = os.path.splitext(in_file)[0] + '.tsv'

    cmd = [
        'msgf_plus',
        'edu.ucsd.msjava.ui.MzIDToTsv',
        '-i', in_file,
        '-o', tmp_file,
        '-showDecoy', 1,
    ]

    cli.execute(*cmd)

    shutil.move(tmp_file, out_file)


def convert_msgf_to_final(in_file, out_file):
    df = load_msgf_df(in_file)

    df.to_csv(out_file, compression='gzip', index=False, sep='\t')


def load_msgf_df(file_name):
    df = pd.read_csv(file_name, sep='\t')

    df = df.rename(columns={
        'SpecID': 'scan_id',
        'ScanNum': 'scan_num',
        'FragMethod': 'frag_method',
        'Precursor': 'precursor_mass',
        'IsotopeError': 'isotope_error',
        'PrecursorError(ppm)': 'precursor_error_ppm',
        'DeNovoScore': 'denovo_score',
        'MSGFScore': 'msgf_score',
        'SpecEValue': 'spec_e_value',
        'EValue': 'e_value',
        'QValue': 'q_value',
        'PepQValue': 'pep_q_value'
    })

    df = df.rename(columns=lambda x: x.lower())

    is_decoy = df['protein'].str.startswith('XXX_').astype(int)

    df['decoy_fdr'] = is_decoy.cumsum() / is_decoy.sum()

    return df


def merge_results(in_files, out_file):
    data = []

    for file_name in soil.utils.workflow.flatten_input(in_files):
        data.append(pd.read_csv(file_name, sep='\t'))

    data = pd.concat(data)

    data = data.sort_values(by='SpecEValue')

    data = data.drop('#SpecFile', axis=1)

    data.to_csv(out_file, index=False, sep='\t')


def run_search_sentinel(sentinel_file, *args, **kwargs):
    db_file = sentinel_file.replace('.sentinel', '.fasta')

    run_search(db_file, *args, **kwargs)


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
        max_mods=1,
        num_tolerable_termini=2,
        num_threads=1,
        precursor_mass_tolerance='20ppm',
        variable_mods=None):

    os.environ['MALLOC_ARENA_MAX'] = '2'

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    mod_file = os.path.join(tmp_dir, 'mods.txt')

    tmp_file = os.path.join(tmp_dir, 'search.mzid')

    write_mods_files(mod_file, max_mods, fixed_mods, variable_mods)

    cmd = [
        'msgf_plus',
        '-Xmx6G',
        '-d', db_file,
        '-s', in_file,
        '-o', tmp_file,
        '-e', enzyme,
        '-addFeatures', int(add_features),
        '-inst', instrument,
        '-m', fragment_method,
        '-mod', mod_file,
        '-ntt', num_tolerable_termini,
        '-t', precursor_mass_tolerance,
        '-tda', int(add_decoys),
        '-thread', num_threads,
        '-ti', ','.join([str(x) for x in isotope_range]),
    ]

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

    if '-term' in aa.lower():
        aa = '*'

    mods_file = soil.utils.package_data.load_data_file('data/protein_mods.tsv')

    mods_df = pd.read_csv(mods_file, sep='\t')

    mods_df = mods_df.set_index('modification')

    row = mods_df.loc[mod]

    weight = row['weight']

    term_spec = row['terminus_specificity']

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
