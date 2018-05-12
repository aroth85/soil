import pipes
import pypeliner.commandline as cli
import subprocess


def run_vardict_paired(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        region,
        out_file,
        min_allele_frequency=0.01,
        remove_duplicate_reads=True):

    cmd = [
        'vardict-java',
        '-b', pipes.quote('{0}|{1}'.format(tumour_bam_file, normal_bam_file)),
        '-f', min_allele_frequency,
        '-G', ref_genome_fasta_file,
        '-R', region,
        '-th', 1
    ]

    if remove_duplicate_reads:
        cmd.append('-t')

    cmd.extend(['>', out_file])

    cmd_str = ' '.join([str(x) for x in cmd])

    subprocess.check_call(cmd_str, shell=True)


def run_test_somatic(in_file, out_file):
    cmd = [
        'cat', in_file, '|', 'testsomatic.R', '>', out_file
    ]

    cli.execute(*cmd)


def run_build_paired_vcf(in_file, out_file):
    cmd = [
        'cat', in_file, '|', 'var2vcf_paired.pl', '>', out_file
    ]

    cli.execute(*cmd)
