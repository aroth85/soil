import os


class SoilRefDataPaths(object):
    """ Simple class to handle access to reference data paths.

    This class should be used to get file paths when using the API.
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    @property
    def cosmic_vcf_file(self):
        """ Path to bgzip compressed tabix indexed COSMIC VCF.

        This file includes the coding and non-coding variants.
        """
        return os.path.join(self.base_dir, 'dbsnp.vcf.gz')

    @property
    def dbsnp_vcf_file(self):
        """ Path to bgzip compressed tabix indexed dbSNP VCF.
        """
        return os.path.join(self.base_dir, 'dbsnp.vcf.gz')

    @property
    def gene_annotations_gtf_file(self):
        """ Path of gene annotation file in GTF format.
        """
        return os.path.join(self.base_dir, 'gene_annotations.gtf')

    @property
    def genome_fasta_file(self):
        """ Path of reference genome FASTA file.
        """
        return os.path.join(self.base_dir, 'genome.fa')

    @property
    def proteome_fasta_file(self):
        """ Path of reference proteome FASTA file.
        """
        return os.path(self.base_dir, 'proteome.fa')

    @property
    def bwa_genome_fasta_file(self):
        """ Path of reference genome FASTA file with BWA index in same directory.

        This is a link to the reference as genome_fasta_file but will be in a separate directory with a BWA index.
        """
        return os.path.join(self.base_dir, 'bwa', 'genome.fa')

    @property
    def star_genome_fasta_file(self):
        """ Path of reference genome FASTA file with STAR index in same directory.

        This is a link to the reference as genome_fasta_file but will be in a separate directory with a STAR index.
        """
        return os.path.join(self.base_dir, 'star', 'genome.fa')
