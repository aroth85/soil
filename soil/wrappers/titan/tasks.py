import glob
import itertools
import numpy as np
import os
import pandas as pd
import pkg_resources
import pypeliner.commandline as cli
import pysam
import shutil
import tarfile
import vcf

import soil.utils.workflow


def split_vcf(in_file, out_file_callback, tmp_dir, split_size=int(1e5)):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_file = os.path.join(tmp_dir, os.path.basename(in_file))

    os.link(in_file, tmp_file)

    cmd = [
        'SnpSift',
        'split',
        '-l', split_size,
        tmp_file,
    ]

    cli.execute(*cmd)

    prefix = tmp_file.split('.')[0]

    os.unlink(tmp_file)

    tmp_split_files = sorted(glob.glob(prefix + '*.vcf'))

    for idx, tmp_out_file in enumerate(tmp_split_files):
        out_file = out_file_callback[idx]

        shutil.move(tmp_out_file, out_file)

    shutil.rmtree(tmp_dir)


def get_snv_allele_counts_for_vcf_targets(
        bam_file,
        vcf_file,
        out_file,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0):

    nucleotides = ('A', 'C', 'G', 'T')

    bam = pysam.AlignmentFile(bam_file, 'rb')

    vcf_reader = vcf.Reader(filename=vcf_file)

    data = []

    for record in vcf_reader:

        counts = bam.count_coverage(
            contig=record.CHROM,
            start=record.POS - 1,
            stop=record.POS,
            quality_threshold=min_bqual,
            read_callback='all'
        )

        counts = np.array(counts, dtype=int).T[0]

        counts = dict(zip(nucleotides, counts))

        ref_base = str(record.REF).upper()

        if len(record.ALT) > 1:
            continue

        alt_base = str(record.ALT[0]).upper()

        # Format output record
        out_row = {
            'chrom': record.CHROM,
            'coord': record.POS,
            'ref': ref_base,
            'ref_counts': counts[ref_base],
            'alt': alt_base,
            'alt_counts': counts[alt_base],
        }

        data.append(out_row)

    data = pd.DataFrame(data, columns=['chrom', 'coord', 'ref', 'ref_counts', 'alt', 'alt_counts'])

    data.to_csv(out_file, index=False, sep='\t')


def merge_counts(in_files, out_file):

    in_files = soil.utils.workflow.flatten_input(in_files)

    pd.read_csv(in_files.pop(0), sep='\t').to_csv(out_file, mode='w', index=False, sep='\t')

    for file_name in in_files:
        pd.read_csv(file_name, sep='\t').to_csv(out_file, mode='a', header=False, index=False, sep='\t')


def build_coverage_file(normal_wig, tumour_wig, gc_wig, mappability_wig, out_file, target_file=None):

    script = pkg_resources.resource_filename('soil', 'scripts/build_titan_coverage.py')

    cmd = [
        'python',
        script,
        '-n', normal_wig,
        '-t', tumour_wig,
        '-g', gc_wig,
        '-m', mappability_wig,
        '-o', out_file,
    ]

    if target_file is not None:
        cmd.extend(['--target-bed-file', target_file])

    cli.execute(*cmd)


def create_intialization_parameters():
    """ Initialize parameter sweep
    """

    normal_contam = [0.25, 0.5, 0.75]

    num_clusters = [1, 2, 3, 4]

    ploidy = [2, 3, 4]

    init_param_values = itertools.product(
        normal_contam,
        num_clusters,
        ploidy,
    )

    init_param_cols = [
        'normal_contam',
        'num_clusters',
        'ploidy',
    ]

    init_params = {}

    for idx, params in enumerate(init_param_values):
        init_params[idx] = dict(zip(init_param_cols, params))

    return init_params


def run_titan(
        coverage_file,
        snp_file,
        init_params,
        out_file,
        tmp_dir,
        estimate_normal_contam=True,
        estimate_ploidy=True,
        is_exome=False,
        sample='Tumour',
        threads=1):

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    script = pkg_resources.resource_filename('soil', 'scripts/run_titan.R')

    cmd = [
        script,
        '--id', sample,
        '--cnFile', coverage_file,
        '--hetFile', snp_file,

        '--outDir', tmp_dir,

        '--normal_0', init_params['normal_contam'],
        '--numClusters', init_params['num_clusters'],
        '--ploidy_0', init_params['ploidy'],

        '--maxCN', 8,
        '--numCores', threads,
    ]

    if estimate_normal_contam:
        cmd.extend(['--estimateNormal', 'map'])
    else:
        cmd.extend(['--estimateNormal', 'fixed'])

    if estimate_ploidy:
        cmd.extend(['--estimatePloidy', 'TRUE'])
    else:
        cmd.extend(['--estimatePloidy', 'FALSE'])

    if is_exome:
        cmd.extend(['--alphaK', 2500])
        cmd.extend(['--alphaKHigh', 2500])

    cli.execute(*cmd)

    with tarfile.open(out_file, 'w:gz') as tar:
        tar.add(tmp_dir, arcname='')

    shutil.rmtree(tmp_dir)


def build_run_stats_file(in_files, init_params, out_file):

    init_params = pd.DataFrame.from_dict(init_params, orient='index')

    init_params = init_params.rename(columns={'normal_contam': 'normal_contam_init', 'ploidy': 'ploidy_init'})

    for init_param_idx in init_params.index:
        titan_params = _read_titan_params(in_files[init_param_idx])

        init_params.loc[init_param_idx, 'log_likelihood'] = titan_params['Log likelihood'][0]

        init_params.loc[init_param_idx, 'normal_contam_est'] = titan_params['Normal contamination estimate'][0]

        init_params.loc[init_param_idx, 'ploidy_est'] = titan_params['Average tumour ploidy estimate'][0]

        init_params.loc[init_param_idx, 's_dbw_validity_index'] = titan_params['S_Dbw validity index (Both)'][0]

    init_params.index.name = 'run_idx'

    init_params = init_params.sort_values(by='s_dbw_validity_index')

    init_params.to_csv(out_file, sep='\t')


def _read_titan_params(tar_file):
    archive = tarfile.open(tar_file, 'r:gz')

    for name in archive.getnames():
        if name.endswith('params.txt'):
            params_name = name

    params = dict()

    params_file = archive.extractfile(params_name)

    for line in params_file:
        key, value = line.split(':')
        params[key] = np.array(value.split()).astype(float)

    archive.close()

    return params


def build_final_results_file(counts_file, coverage_file, run_files, stats_file, out_file, tmp_dir):

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    input_dir = os.path.join(tmp_dir, 'input')

    os.makedirs(input_dir)

    shutil.copyfile(stats_file, os.path.join(tmp_dir, 'stats.tsv'))

    shutil.copyfile(counts_file, os.path.join(input_dir, 'counts.tsv'))

    shutil.copyfile(coverage_file, os.path.join(input_dir, 'coverage.tsv'))

    for run_idx in run_files:
        run_tmp_dir = os.path.join(tmp_dir, 'runs', str(run_idx))

        os.makedirs(run_tmp_dir)

        tar_file = tarfile.open(run_files[run_idx], 'r:gz')

        tar_file.extractall(run_tmp_dir)

        _rename_titan_files(run_tmp_dir)

    df = pd.read_csv(stats_file, sep='\t')

    df = df.sort_values(by='s_dbw_validity_index')

    best_idx = df.iloc[0]['run_idx']

    shutil.copytree(os.path.join(tmp_dir, 'runs', str(int(best_idx))), os.path.join(tmp_dir, 'runs', 'selected'))

    with tarfile.open(out_file, 'w:gz') as tar:
        tar.add(tmp_dir, arcname='')

    shutil.rmtree(tmp_dir)


def _rename_titan_files(run_dir):
    plot_dir = os.path.join(run_dir, 'plots')

    prefix = os.path.basename(os.path.commonprefix(os.listdir(run_dir)))

    shutil.move(os.path.join(run_dir, prefix), plot_dir)

    shutil.move(os.path.join(run_dir, prefix + '.params.txt'), os.path.join(run_dir, 'params.txt'))

    shutil.move(os.path.join(run_dir, prefix + '.seg'), os.path.join(run_dir, 'igv.seg'))

    shutil.move(os.path.join(run_dir, prefix + '.segs.txt'), os.path.join(run_dir, 'segs.txt'))

    shutil.move(os.path.join(run_dir, prefix + '.titan.txt'), os.path.join(run_dir, 'pos.txt'))

    shutil.move(os.path.join(run_dir, prefix + '.RData'), os.path.join(run_dir, 'workspace.RData'))

    for file_name in glob.glob(os.path.join(plot_dir, '*')):
        shutil.move(file_name, os.path.join(plot_dir, os.path.basename(file_name).replace(prefix + '_', '')))
