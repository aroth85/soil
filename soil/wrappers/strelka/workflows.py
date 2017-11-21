import pandas as pd
import pypeliner
import pypeliner.managed as mgd
import pysam

import soil.utils.genome
import soil.utils.workflow
import soil.wrappers.samtools.tasks
import soil.wrappers.strelka.tasks

low_mem_ctx = {'mem': 2, 'mem_retry_factor': 2, 'num_retry': 3}
med_mem_ctx = {'mem': 4, 'mem_retry_factor': 2, 'num_retry': 3}


def create_somatic_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        chromosomes=None,
        is_exome=False,
        split_size=int(1e7)):

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools', 'strelka'])

    workflow = pypeliner.workflow.Workflow(default_ctx=med_mem_ctx, default_sandbox=sandbox)

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
        kwargs={
            'is_exome': is_exome,
        }
    )

    workflow.transform(
        name='merge_indels',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('indels.vcf', 'regions'),
            mgd.TempOutputFile('indels.vcf.gz'),
        ),
    )

    workflow.transform(
        name='merge_snvs',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('snvs.vcf', 'regions'),
            mgd.TempOutputFile('snvs.vcf.gz'),
        ),
    )

    workflow.transform(
        name='merge_all',
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
        ctx=low_mem_ctx,
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
        ctx=low_mem_ctx,
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
