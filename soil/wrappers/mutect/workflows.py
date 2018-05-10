import pysam
import pypeliner
import pypeliner.managed as mgd

import soil.utils.genome
import soil.utils.workflow
import soil.wrappers.samtools.tasks

import tasks

low_mem_ctx = {'mem': 2, 'mem_retry_factor': 2, 'num_retry': 3}
med_mem_ctx = {'mem': 4, 'mem_retry_factor': 2, 'num_retry': 3}


def create_mutect_paired_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        chromosomes=None,
        normal_name='normal',
        split_size=int(1e7),
        tumour_name='tumour'):

    normal_name = get_sample(normal_bam_file, normal_name)

    tumour_name = get_sample(tumour_bam_file, tumour_name)

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'gatk', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_ctx=low_mem_ctx, default_sandbox=sandbox)

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('config', 'regions'),
        value=soil.utils.genome.get_bam_regions(normal_bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.transform(
        name='run_mutect',
        axes=('regions',),
        ctx=med_mem_ctx,
        func=tasks.run_mutect,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempInputObj('config', 'regions'),
            mgd.TempOutputFile('region.vcf.gz', 'regions')
        ),
        kwargs={
            'normal_name': normal_name,
            'tumour_name': tumour_name
        }
    )

    workflow.transform(
        name='concatenate_vcfs',
        func=soil.wrappers.samtools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('region.vcf.gz', 'regions'),
            mgd.OutputFile(out_file),
        ),
    )

    return workflow


def get_sample(file_name, orig_name):
    bam = pysam.AlignmentFile(file_name)

    try:
        sample = bam.header['RG'][0]['SM']

    except KeyError:
        sample = orig_name

    bam.close()

    return sample
