class Package(object):

    def __init__(self, name, version, channel='bioconda'):
        self.name = name

        self.version = version

        self.channel = channel

    @property
    def conda_str(self):
        return '{0} =={1}'.format(self.name, self.version)


packages = {
    'bwa': Package('bwa', '0.7.16'),
    'bcftools': Package('bcftools', '1.6'),
    'opossum': Package('opossum', '0.2.0de259b45a35cd7d4c01dbb40d48d16c9118d76c4', channel='aroth85'),
    'optitype': Package('optitype', '1.2.1'),
    'mixcr': Package('mixcr', '2.1.3'),
    'msgf_plus': Package('msgf_plus', '2017.07.21'),
    'percolator': Package('percolator', '3.1'),
    'platypus': Package('platypus-variant', '0.8.1.1'),
    'proteowizard': Package('proteowizard', '3_0_9992'),
    'sambamba': Package('sambamba', '0.6.6'),
    'samtools': Package('samtools', '1.6'),
    'star': Package('star', '2.5.3a'),
    'strelka': Package('strelka', '2.8.4'),
    'stringtie': Package('stringtie', '1.3.3'),
    'transdecoder': Package('transdecoder', '3.0.1'),
    'varscan': Package('varscan', '2.4.3'),
}
