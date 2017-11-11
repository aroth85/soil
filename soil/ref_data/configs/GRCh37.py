def get_config():
    cosmic = {
        'host': 'sftp-cancer.sanger.ac.uk',
        'cosmic_version': 'v83'
    }

    dbsnp_url = 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz'

    gencode_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping'

    ref_genome_fasta_url = _get_url(gencode_url, 'GRCh37.primary_assembly.genome.fa.gz')

    ref_proteome_fasta_url = _get_url(gencode_url, 'gencode.v27lift37.pc_translations.fa.gz')

    ref_gene_annotations_gtf_url = _get_url(gencode_url, 'gencode.v27lift37.annotation.gtf.gz')

    config = locals()

    return config


def _get_url(base_url, file_name):
    return '/'.join([base_url, file_name])
