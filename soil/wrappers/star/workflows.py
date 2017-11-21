import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.sambamba.tasks

import tasks


def create_multiple_lane_align_workflow(
        fastq_files_1,
        fastq_files_2,
        ref_genome_dir,
        out_bam_file,
        add_xs_tag=False,
        align_threads=1,
        merge_threads=1,
        read_group_info=None,
        sort_threads=1):

    if read_group_info is None:
        read_group_info = {}

        for key in fastq_files_1:
            read_group_info[key] = None

    sandbox = soil.utils.workflow.get_sandbox(['sambamba'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(obj=mgd.OutputChunks('lane'), value=fastq_files_1.keys())

    workflow.subworkflow(
        name='align',
        axes=('lane',),
        func=create_align_workflow,
        args=(
            mgd.InputFile('R1.fq.gz', 'lane', fnames=fastq_files_1),
            mgd.InputFile('R2.fq.gz', 'lane', fnames=fastq_files_2),
            ref_genome_dir,
            mgd.TempOutputFile('lane.bam', 'lane'),
        ),
        kwargs={
            add_xs_tag: add_xs_tag,
            align_threads: align_threads,
            read_group_info: read_group_info,
            sort_threads: sort_threads,
        }
    )

    workflow.transform(
        name='markdups_and_merge',
        axes=(),
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': merge_threads},
        func=soil.wrappers.sambamba.tasks.markdups,
        args=(
            mgd.TempInputFile('lane.bam', 'lane'),
            mgd.OutputFile(out_bam_file),
            mgd.TempSpace('markdup_tmp'),
        ),
        kwargs={
            'threads': merge_threads,
        }
    )

    return workflow


def create_align_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_dir,
        out_bam_file,
        add_xs_tag=False,
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
            'add_xs_tag': add_xs_tag,
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
