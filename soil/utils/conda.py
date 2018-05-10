class Package(object):

    def __init__(self, name, version=None, channel='bioconda'):
        self.name = name

        self.version = version

        self.channel = channel

    @property
    def conda_str(self):
        if self.version is None:
            pkg_str = self.name

        else:
            pkg_str = '{0} =={1}'.format(self.name, self.version)

        return pkg_str


packages = {
    # B
    'bwa': Package('bwa', '0.7.16'),
    'bcftools': Package('bcftools', '1.6'),
    # E
    'eagle': Package('eagle-phase', '2.3.5', channel='soil'),
    # G
    'gatk': Package('gatk4', '4.0.4.0'),
    # H
    'hmmcopy': Package('bioconductor-hmmcopy', '1.20.0'),
    'hmmcopy_utils': Package('hmmcopy', '0.1.1'),
    # K
    'kallisto': Package('kallisto', '0.43.1'),
    # M
    'mixcr': Package('mixcr', '2.1.3'),
    'msgf_plus': Package('msgf_plus', '2017.07.21'),
    # O
    'opossum': Package('opossum', '0.2.0de259b45a35cd7d4c01dbb40d48d16c9118d76c4', channel='aroth85'),
    'optitype': Package('optitype', '1.2.1'),
    # P
    'pyensembl': Package('pyensembl', '1.1.0', channel='pip'),
    'percolator': Package('percolator', '3.1'),
    'platypus': Package('platypus-variant', '0.8.1.1'),
    'proteowizard': Package('proteowizard', '3_0_9992'),
    # R
    'razers3': Package('razers3', '3.5.0'),
    # S
    'sambamba': Package('sambamba', '0.6.6'),
    'samtools': Package('samtools', '1.6'),
    'snpeff': Package('snpeff', '4.3.1r'),
    'snpsift': Package('snpsift', '4.3.1r'),
    'star': Package('star', '2.5.3a'),
    'strelka': Package('strelka', '2.8.4'),
    'stringtie': Package('stringtie', '1.3.3'),
    'subread': Package('subread', '1.6.0'),
    # T
    'topiary': Package('topiary', '2.1.0', channel='pip'),
    'titan': Package('bioconductor-titancna', '1.16.0'),
    'transdecoder': Package('transdecoder', '3.0.1'),
    # U
    'ucsc-bedgraphtobigwig': Package('ucsc-bedgraphtobigwig', '357'),
    # V
    'varscan': Package('varscan', '2.4.3'),
}
