from pypeliner.sandbox import CondaSandbox

import pandas as pd
import pypeliner
import pypeliner.managed as mgd
import pysam

import soil.wrappers.samtools.tasks
import soil.wrappers.strelka.tasks
import soil.utils.genome


def create_somatic_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        chromosomes=None,
        split_size=int(1e7)):

    sandbox = CondaSandbox(
        channels=('bioconda',),
        packages=('bcftools ==1.6', 'samtools ==1.6', 'strelka ==2.8.4'),
    )

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(
        obj=mgd.TempOutputObj('config', 'regions'),
        value=soil.utils.genome.get_bam_regions(normal_bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.setobj(
        obj=mgd.TempOutputObj('chrom_names', 'chrom_axis'),
        value=get_chromosomes(normal_bam_file, chromosomes=chromosomes)
    )

    workflow.transform(
        name='count_fasta_bases',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=soil.wrappers.strelka.tasks.count_fasta_bases,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('ref_base_counts.tsv'),
        ),
    )

    workflow.transform(
        name='get_genome_size',
        ctx={'local': True},
        func=get_known_genome_size,
        ret=mgd.TempOutputObj('genome_size'),
        args=(
            mgd.InputFile(tumour_bam_file),
            mgd.TempInputFile('ref_base_counts.tsv'),
            chromosomes,
        ),
        sandbox=None,
    )

    workflow.transform(
        name='get_chromosome_depths',
        axes=('chrom_axis',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 8},
        func=soil.wrappers.strelka.tasks.get_chromosome_depth,
        args=(
            mgd.TempInputObj('chrom_names', 'chrom_axis'),
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('chrom_depth.txt', 'chrom_axis'),
        ),
    )

    workflow.transform(
        name='merge_chromosome_depths',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 8},
        func=soil.wrappers.strelka.tasks.merge_chromosome_depth,
        args=(
            mgd.TempInputFile('chrom_depth.txt', 'chrom_axis'),
            mgd.TempOutputFile('chrom_depth_merged.txt'),
        ),
        sandbox=None,
    )

    workflow.transform(
        name='call_genome_segment',
        axes=('regions',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 8},
        func=soil.wrappers.strelka.tasks.call_genome_segment,
        args=(
            mgd.TempInputFile('chrom_depth_merged.txt'),
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('indels.vcf', 'regions'),
            mgd.TempOutputFile('snvs.vcf', 'regions'),
            mgd.TempSpace('call_genome_segment_tmp', 'regions'),
            mgd.TempInputObj('config', 'regions'),
            mgd.TempInputObj('genome_size'),
        ),
    )

    workflow.transform(
        name='merge_indels',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('indels.vcf', 'regions'),
            mgd.TempOutputFile('indels.vcf.gz'),
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('snvs.vcf', 'regions'),
            mgd.TempOutputFile('snvs.vcf.gz'),
        ),
    )

    workflow.transform(
        name='merge_all',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            [mgd.TempInputFile('indels.vcf.gz'), mgd.TempInputFile('snvs.vcf.gz')],
            mgd.TempOutputFile('merged.vcf.gz'),
        ),
        kwargs={
            'allow_overlap': True,
        },
    )

    workflow.commandline(
        name='filter_vcf',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        args=(
            'bcftools',
            'view',
            '-O', 'z',
            '-f', '.,PASS',
            '-o', mgd.OutputFile(out_file),
            mgd.TempInputFile('merged.vcf.gz'),
        )
    )

    workflow.transform(
        name='index_vcf',
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=soil.wrappers.samtools.tasks.index_vcf,
        args=(
            mgd.InputFile(out_file),
            mgd.OutputFile(out_file + '.tbi'),
        )
    )

    return workflow


def get_chromosomes(bam_file, chromosomes=None):
    chroms = {}

    for idx, chrom in enumerate(_get_chromosomes(bam_file, chromosomes=chromosomes)):
        chroms[idx] = chrom

    return chroms


def get_known_genome_size(bam_file, size_file, chromosomes):
    chromosomes = _get_chromosomes(bam_file, chromosomes)

    sizes = pd.read_csv(
        size_file,
        converters={'chrom': str},
        header=None,
        names=['path', 'chrom', 'known_size', 'size'],
        sep='\t'
    )

    sizes = sizes[sizes['chrom'].isin(chromosomes)]

    return sizes['known_size'].sum()


def _get_chromosomes(bam_file, chromosomes=None):
    bam = pysam.Samfile(bam_file, 'rb')

    if chromosomes is None:
        chromosomes = bam.references

    else:
        chromosomes = chromosomes

    return [str(x) for x in chromosomes]
