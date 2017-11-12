import os
import pypeliner.commandline as cli
import shutil

import soil.utils.file_system


def gtf_to_fasta(
        in_gtf_file,
        ref_genome_fasta_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_genome_gff_file,
        tmp_dir):
    """ Convert a GTF file from cufflinks or stringtie to a FASTA file
    """
    os.chdir(tmp_dir)

    opt_dir = os.path.join(os.environ['CONDA_PREFIX'], 'opt')

    gtf_to_fasta_script = soil.utils.file_system.find('cufflinks_gtf_genome_to_cdna_fasta.pl', opt_dir)

    gtf_to_gff_script = soil.utils.file_system.find('cufflinks_gtf_to_alignment_gff3.pl', opt_dir)

    cdna_to_genome_script = soil.utils.file_system.find('cdna_alignment_orf_to_genome_orf.pl', opt_dir)

    tmp_cdna_fasta_file = os.path.join(tmp_dir, 'cdna.fasta')

    tmp_gff_file = os.path.join(tmp_dir, 'raw.gff3')

    cmd = ['perl', gtf_to_fasta_script, in_gtf_file, ref_genome_fasta_file, '>', tmp_cdna_fasta_file]

    cli.execute(*cmd)

    cmd = ['perl', gtf_to_gff_script, in_gtf_file, '>', tmp_gff_file]

    cmd = ['TransDecoder.LongOrfs', '-t', tmp_cdna_fasta_file]

    cli.execute(*cmd)

    cmd = ['TransDecoder.Predict', '-t', tmp_cdna_fasta_file]

    cli.execute(*cmd)

    cmd = [
        'perl',
        cdna_to_genome_script,
        tmp_cdna_fasta_file + '.transdecoder.gff3',
        tmp_gff_file,
        tmp_cdna_fasta_file,
        '>',
        out_genome_gff_file
    ]

    cli.execute(*cmd)

    shutil.move(tmp_cdna_fasta_file + '.transdecoder.cds', out_cds_fasta_file)

    shutil.move(tmp_cdna_fasta_file + '.transdecoder.pep', out_protein_fasta_file)
