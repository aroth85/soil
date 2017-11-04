from pypeliner.sandbox import CondaSandbox

import pypeliner
import pypeliner.managed as mgd

import pot.wrappers.bwa
import pot.wrappers.sambamba


def create_single_lane_alignment_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_fasta_file,
        out_bam_file,
        bwa_threads=1,
        read_group_info=None,
        sambamba_threads=1):

    sandbox = CondaSandbox(
        channels=('bioconda',),
        packages=('bwa ==0.7.16', 'samtools ==1.6', 'sambamba ==0.6.6'),
    )

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='bwa_mem_paired_end',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': bwa_threads},
        func=pot.wrappers.bwa.mem_paired_end,
        args=(
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('aligned.bam'),
        ),
        kwargs={
            'read_group_info': read_group_info,
            'threads': bwa_threads,
        }
    )

    workflow.transform(
        name='sort',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': sambamba_threads},
        func=pot.wrappers.sambamba.sort,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.OutputFile(out_bam_file),
            mgd.TempSpace('sort_tmp'),
        ),
        kwargs={'threads': sambamba_threads}
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
