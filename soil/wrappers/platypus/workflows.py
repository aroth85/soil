import pypeliner
import pypeliner.managed as mgd

import soil.utils.genome
import soil.utils.workflow
import soil.wrappers.samtools.tasks

import tasks

low_mem_ctx = {'mem': 2, 'mem_retry_factor': 2, 'num_retry': 3}
med_mem_ctx = {'mem': 4, 'mem_retry_factor': 2, 'num_retry': 3}


def create_single_sample_workflow(
        bam_file,
        ref_genome_fasta_file,
        out_file, chromosomes='default',
        split_size=int(1e7)):

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools', 'platypus'])

    workflow = pypeliner.workflow.Workflow(default_ctx=low_mem_ctx, default_sandbox=sandbox)

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('config', 'regions'),
        value=soil.utils.genome.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.transform(
        name='run_platypus',
        axes=('regions',),
        ctx=med_mem_ctx,
        func=tasks.call_variants,
        args=(
            mgd.InputFile(bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('log.txt', 'regions'),
            mgd.TempOutputFile('region.vcf', 'regions'),
            mgd.TempInputObj('config', 'regions'),
        ),
    )

    workflow.transform(
        name='compress',
        axes=('regions',),
        func=soil.wrappers.samtools.tasks.compress_vcf,
        args=(
            mgd.TempInputFile('region.vcf', 'regions'),
            mgd.TempOutputFile('region.vcf.gz', 'regions'),
        ),
    )

    workflow.transform(
        name='concatenate_vcfs',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('region.vcf.gz', 'regions'),
            mgd.TempOutputFile('merged.vcf.gz'),
        ),
    )

    workflow.commandline(
        name='filter_vcf',
        args=(
            'bcftools',
            'view',
            '-O', 'z',
            '-f', '.,PASS',
            '-o', mgd.OutputFile(out_file),
            mgd.TempInputFile('merged.vcf.gz'),
        )
    )
    return workflow


def create_rna_single_sample_workflow(
        bam_file,
        ref_genome_fasta_file,
        out_file,
        chromosomes='default',
        split_size=int(1e7)):

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools', 'opossum', 'platypus'])

    workflow = pypeliner.workflow.Workflow(default_ctx=low_mem_ctx, default_sandbox=sandbox)

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('config', 'regions'),
        value=soil.utils.genome.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.commandline(
        name='split_bam',
        axes=('regions',),
        args=(
            'samtools',
            'view',
            '-b',
            mgd.InputFile(bam_file),
            mgd.TempInputObj('config', 'regions'),
            '>',
            mgd.TempOutputFile('raw.bam', 'regions'),
        )
    )

    workflow.transform(
        name='run_opossum',
        axes=('regions',),
        ctx=med_mem_ctx,
        func=tasks.fix_rna_bam,
        args=(
            mgd.TempInputFile('raw.bam', 'regions'),
            mgd.TempOutputFile('fixed.bam', 'regions'),
        )
    )

    workflow.transform(
        name='run_platypus',
        axes=('regions',),
        ctx=med_mem_ctx,
        func=tasks.call_variants,
        args=(
            mgd.TempInputFile('fixed.bam', 'regions'),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('log.txt', 'regions'),
            mgd.TempOutputFile('region.vcf', 'regions'),
            mgd.TempInputObj('config', 'regions'),
        ),
        kwargs={
            'rna': True
        }
    )

    workflow.transform(
        name='compress',
        axes=('regions',),
        func=soil.wrappers.samtools.tasks.compress_vcf,
        args=(
            mgd.TempInputFile('region.vcf', 'regions'),
            mgd.TempOutputFile('region.vcf.gz', 'regions'),
        ),
    )

    workflow.transform(
        name='concatenate_vcfs',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('region.vcf.gz', 'regions'),
            mgd.TempOutputFile('merged.vcf.gz'),
        ),
    )

    workflow.commandline(
        name='filter_vcf',
        args=(
            'bcftools',
            'view',
            '-O', 'z',
            '-f', '.,PASS',
            '-o', mgd.OutputFile(out_file),
            mgd.TempInputFile('merged.vcf.gz'),
        )
    )
    return workflow
