import pipes
import pypeliner
import pypeliner.managed as mgd

import soil.utils.genome
import soil.utils.workflow
import soil.wrappers.samtools.tasks

import tasks

low_mem_ctx = {'mem': 2, 'mem_retry_factor': 2, 'num_retry': 3}
med_mem_ctx = {'mem': 8, 'mem_retry_factor': 2, 'num_retry': 3}


def create_vardict_paired_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        chromosomes=None,
        split_size=int(5e6)):

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools', 'vardict', 'vardict-java'])

    workflow = pypeliner.workflow.Workflow(default_ctx=low_mem_ctx, default_sandbox=sandbox)

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('config', 'regions'),
        value=soil.utils.genome.get_bam_regions(normal_bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.transform(
        name='run_vardict',
        axes=('regions',),
        ctx=med_mem_ctx,
        func=tasks.run_vardict_paired,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempInputObj('config', 'regions'),
            mgd.TempOutputFile('call.tsv', 'regions')
        )
    )

    workflow.transform(
        name='test_somatic',
        axes=('regions',),
        func=tasks.run_test_somatic,
        args=(
            mgd.TempInputFile('call.tsv', 'regions'),
            mgd.TempOutputFile('somatic.tsv', 'regions')
        )
    )

    workflow.transform(
        name='write_vcf',
        axes=('regions',),
        func=tasks.run_build_paired_vcf,
        args=(
            mgd.TempInputFile('somatic.tsv', 'regions'),
            mgd.TempOutputFile('region.vcf', 'regions')
        )
    )

    workflow.commandline(
        name='compress_vcf',
        axes=('regions',),
        args=(
            'bcftools', 'view',
            '-O', 'z',
            '-o', mgd.TempOutputFile('region.vcf.gz', 'regions'),
            mgd.TempInputFile('region.vcf', 'regions')
        )
    )

    workflow.transform(
        name='concatenate_vcfs',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('region.vcf.gz', 'regions'),
            mgd.TempOutputFile('merged.vcf.gz'),
        )
    )

    workflow.commandline(
        name='filter_vcf',
        args=(
            'bcftools',
            'view',
            '-O', 'z',
            '-f', '.,PASS',
            '-o', mgd.TempOutputFile('filtered.vcf.gz'),
            mgd.TempInputFile('merged.vcf.gz'),
        )
    )

    workflow.commandline(
        name='filter_somatics',
        args=(
            'bcftools',
            'filter',
            '-i', pipes.quote('INFO/STATUS[0]="StrongSomatic"'),
            '-O', 'z',
            '-o', mgd.OutputFile(out_file),
            mgd.TempInputFile('filtered.vcf.gz')
        )
    )

    return workflow
