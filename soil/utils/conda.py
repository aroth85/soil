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
    'platypus': Package('platypus-variant', '0.8.1.1'),
    'sambamba': Package('sambamba', '0.6.6'),
    'samtools': Package('samtools', '1.6'),
    'strelka': Package('strelka', '2.8.4'),
    'varscan': Package('varscan', '2.4.3'),
}