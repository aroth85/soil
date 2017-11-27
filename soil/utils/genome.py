from collections import OrderedDict

import pysam


def get_regions(chromosome_lengths, split_size):
    """ Split up chromosomes into regions with a fixed size. Useful for parallelising tasks across a genome.

    :param chromosome_lengths: Dictionary with chromosomes as keys and lenghts as values
    :param split_size: Maximum lenght of split interval.
    :returns: A dictionary with keys being the numeric id of the region and values being the samtools style region
        string.

    >>> get_regions({'1': 1000, '2': 100}, 500)
    {0: '1:1-500', 1: '1:501-1000', 2: '2:1-100'}
    """
    if split_size is None:
        return dict(enumerate(chromosome_lengths.keys()))

    regions = {}

    region_index = 0

    for chrom, length in chromosome_lengths.iteritems():
        lside_interval = range(1, length + 1, split_size)

        rside_interval = range(split_size, length + split_size, split_size)

        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)

            regions[region_index] = '{}:{}-{}'.format(chrom, beg, end)

            region_index += 1

    return regions


def get_bam_regions(bam_file, split_size, chromosomes='default'):
    chromosome_lengths = load_bam_chromosome_lengths(bam_file, chromosomes=chromosomes)

    return get_regions(chromosome_lengths, split_size)


def get_default_chromosomes(bam_file):
    bam = pysam.Samfile(bam_file, 'rb')

    defaults = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']

    chromosomes = []

    for chrom in bam.references:
        raw_chrom = chrom.replace('chr', '')

        if raw_chrom in defaults:
            chromosomes.append(chrom)

    return chromosomes


def load_bam_chromosome_lengths(file_name, chromosomes='default'):
    chromosome_lengths = OrderedDict()

    bam = pysam.Samfile(file_name, 'rb')

    if chromosomes == 'all':
        chromosomes = bam.references

    elif chromosomes == 'default':
        chromosomes = get_default_chromosomes(file_name)

    else:
        chromosomes = chromosomes

    for chrom, length in zip(bam.references, bam.lengths):
        if chrom not in chromosomes:
            continue

        chromosome_lengths[str(chrom)] = int(length)

    return chromosome_lengths
