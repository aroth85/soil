import pypeliner.commandline as cli


def run_mutect(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        region,
        out_file,
        normal_name='normal',
        tumour_name='tumour'):

    cmd = [
        'gatk',
        'Mutect2',
        '-R', ref_genome_fasta_file,
        '-I', tumour_bam_file,
        '-tumor', tumour_name,
        '-I', normal_bam_file,
        '-normal', normal_name,
        '-L', region,
        '-O', out_file
    ]

    cli.execute(*cmd)
