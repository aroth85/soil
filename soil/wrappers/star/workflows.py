import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.sambamba.tasks

import tasks


def create_align_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_dir,
        out_bam_file,
        align_threads=1,
        read_group_info=None,
        sort_threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['star', 'samtools', 'sambamba'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='star_align',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': align_threads},
        func=tasks.align,
        args=(
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            ref_genome_dir,
            mgd.TempOutputFile('aligned.bam'),
            mgd.TempSpace('align_tmp'),
        ),
        kwargs={
            'read_group_info': read_group_info,
            'threads': align_threads,
        }
    )

    workflow.transform(
        name='sort',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': sort_threads},
        func=soil.wrappers.sambamba.tasks.sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.OutputFile(out_bam_file),
            mgd.TempSpace('sort_tmp'),
        ),
        kwargs={'threads': sort_threads}
    )

    workflow.commandline(
        name='index',
        args=(
            'samtools',
            'index',
            mgd.InputFile(out_bam_file),
            mgd.OutputFile(out_bam_file + '.bai'),
        )
    )

    return workflow
