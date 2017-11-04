"""
Wrappers for tools included with the Strelka program. Note that to actually use strelka you will need to use one of the
workflows.

https://github.com/Illumina/strelka
"""
import ConfigParser
import os
import pypeliner.commandline as cli
import shutil

import pot.utils.file_system


def call_genome_segment(
        chrom_depth_file,
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_file,
        snv_file,
        tmp_dir,
        region,
        genome_size,
        is_exome=False):

    def load_config(file_name):
        config_parser = ConfigParser.ConfigParser()

        config_parser.optionxform = str

        config_parser.read(file_name)

        return dict(config_parser.items('StrelkaSomatic'))

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_indel_file = os.path.join(tmp_dir, 'indels.vcf')

    tmp_snv_file = os.path.join(tmp_dir, 'snvs.vcf')

    stats_file = os.path.join(tmp_dir, 'stats.txt')

    share_dir = os.path.join(os.environ['CONDA_PREFIX'], 'share')

    config = load_config(pot.utils.file_system.find('configureStrelkaSomaticWorkflow.py.ini', share_dir))

    strelka_exe = pot.utils.file_system.find('strelka2', share_dir)

    cmd = [
        strelka_exe,
        '--normal-align-file', normal_bam_file,
        '--tumor-align-file', tumour_bam_file,

        '--somatic-indel-file', tmp_indel_file,
        '--somatic-snv-file', tmp_snv_file,
        '--stats-file', stats_file,

        # strelkaSharedWorkflow.py
        '--region', region,
        '--ref', ref_genome_fasta_file,
        '-genome-size', genome_size,
        '-max-indel-size', 50,

        # strelkaSomaticWorkflow.py
        '-min-mapping-quality', config['minTier1Mapq'],
        '-min-qscore', 0,
        '-max-window-mismatch', 3, 20,
        '-indel-nonsite-match-prob', 0.5,
        '--somatic-snv-rate', config['ssnvPrior'],
        '--shared-site-error-rate', config['ssnvNoise'],
        '--shared-site-error-strand-bias-fraction', config['ssnvNoiseStrandBiasFrac'],
        '--somatic-indel-rate', config['sindelPrior'],
        '--shared-indel-error-factor', config['sindelNoiseFactor'],
        '--tier2-min-mapping-quality', config['minTier2Mapq'],
        '--tier2-mismatch-density-filter-count', 10,
        '--tier2-indel-nonsite-match-prob', 0.25,
        '--tier2-include-singleton',
        '--tier2-include-anomalous',

        '--strelka-snv-max-filtered-basecall-frac', config['snvMaxFilteredBasecallFrac'],
        '--strelka-snv-max-spanning-deletion-frac', config['snvMaxSpanningDeletionFrac'],
        '--strelka-snv-min-qss-ref', config['ssnvQuality_LowerBound'],

        '--strelka-indel-max-window-filtered-basecall-frac', config['indelMaxWindowFilteredBasecallFrac'],
        '--strelka-indel-min-qsi-ref', config['sindelQuality_LowerBound'],

        '--ssnv-contam-tolerance', config['ssnvContamTolerance'],
        '--indel-contam-tolerance', config['indelContamTolerance'],

        '--somatic-snv-scoring-model-file', pot.utils.file_system.find('somaticSNVScoringModels.json', share_dir),
        '--somatic-indel-scoring-model-file', pot.utils.file_system.find('somaticIndelScoringModels.json', share_dir),
    ]

    if not is_exome:
        cmd.extend([
            '--strelka-chrom-depth-file', chrom_depth_file,
            '--strelka-max-depth-factor', config['depthFilterMultiple'],
        ])

    cli.execute(*cmd)


def count_fasta_bases(ref_genome_fasta_file, out_file):
    share_dir = os.path.join(os.environ['CONDA_PREFIX'], 'share')

    exe = pot.utils.file_system.find('countFastaBases', share_dir)

    cmd = [
        exe,
        ref_genome_fasta_file,
        '>',
        out_file
    ]

    cli.execute(*cmd)


def get_chromosome_depth(chrom, bam_file, out_file):
    share_dir = os.path.join(os.environ['CONDA_PREFIX'], 'share')

    exe = pot.utils.file_system.find('GetChromDepth', share_dir)

    cmd = [
        exe,
        '--align-file', bam_file,
        '--chrom', chrom,
        '--output-file', out_file
    ]

    cli.execute(*cmd)


def merge_chromosome_depth(in_files, out_file):
    with open(out_file, 'w') as out_fh:
        for key in sorted(in_files):
            with open(in_files[key], 'r') as in_fh:
                for line in in_fh:
                    out_fh.write(line)
