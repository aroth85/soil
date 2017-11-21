import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.bwa.tasks
import soil.wrappers.sambamba.tasks


def create_multiple_lane_align_workflow(
        fastq_files_1,
        fastq_files_2,
        ref_genome_file,
        out_bam_file,
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

    workflow.setobj(obj=mgd.TempOutputObj('read_group_info', 'lane'), value=read_group_info)

    workflow.subworkflow(
        name='align',
        axes=('lane',),
        func=create_align_workflow,
        args=(
            mgd.InputFile('R1.fq.gz', 'lane', fnames=fastq_files_1),
            mgd.InputFile('R2.fq.gz', 'lane', fnames=fastq_files_2),
            mgd.InputFile(ref_genome_file),
            mgd.TempOutputFile('lane.bam', 'lane'),
        ),
        kwargs={
            'align_threads': align_threads,
            'read_group_info': mgd.TempInputObj('read_group_info', 'lane'),
            'sort_threads': sort_threads,
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
        ref_genome_fasta_file,
        out_bam_file,
        align_threads=1,
        read_group_info=None,
        sort_threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['bwa', 'samtools', 'sambamba'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='bwa_mem_paired_end',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': align_threads},
        func=soil.wrappers.bwa.tasks.mem_paired_end,
        args=(
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('aligned.bam'),
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
