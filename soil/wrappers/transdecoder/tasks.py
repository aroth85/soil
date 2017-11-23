import os
import pypeliner.commandline as cli
import shutil

import soil.utils.file_system


def convert_gtf_to_cdna_fasta(in_file, ref_genome_fasta_file, out_file):
    opt_dir = os.path.join(os.environ['CONDA_PREFIX'], 'opt')

    gtf_to_fasta_script = soil.utils.file_system.find('cufflinks_gtf_genome_to_cdna_fasta.pl', opt_dir)

    cmd = ['perl', gtf_to_fasta_script, in_file, ref_genome_fasta_file, '>', out_file]

    cli.execute(*cmd)


def convert_gtf_to_gff_file(in_file, out_file):
    opt_dir = os.path.join(os.environ['CONDA_PREFIX'], 'opt')

    gtf_to_gff_script = soil.utils.file_system.find('cufflinks_gtf_to_alignment_gff3.pl', opt_dir)

    cmd = ['perl', gtf_to_gff_script, in_file, '>', out_file]

    cli.execute(*cmd)


def run_transdecoder(in_file, out_cds_fasta_file, out_protein_fasta_file, out_gff_file, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_in_file = os.path.join(tmp_dir, 'cdna.fasta')

    shutil.copyfile(in_file, tmp_in_file)

    old_wd = os.getcwd()

    os.chdir(tmp_dir)

    cmd = ['TransDecoder.LongOrfs', '-t', tmp_in_file]

    cli.execute(*cmd)

    cmd = ['TransDecoder.Predict', '-t', tmp_in_file]

    cli.execute(*cmd)

    os.chdir(old_wd)

    shutil.move(tmp_in_file + '.transdecoder.cds', out_cds_fasta_file)

    shutil.move(tmp_in_file + '.transdecoder.pep', out_protein_fasta_file)

    shutil.move(tmp_in_file + '.transdecoder.gff3', out_gff_file)

    shutil.rmtree(tmp_dir)


def build_alignment_gff(cdna_fasta_file, gff_file, ref_gff_file, out_file):
    opt_dir = os.path.join(os.environ['CONDA_PREFIX'], 'opt')

    cdna_to_genome_script = soil.utils.file_system.find('cdna_alignment_orf_to_genome_orf.pl', opt_dir)

    cmd = [
        'perl',
        cdna_to_genome_script,
        gff_file,
        ref_gff_file,
        cdna_fasta_file,
        '>',
        out_file
    ]

    cli.execute(*cmd)
