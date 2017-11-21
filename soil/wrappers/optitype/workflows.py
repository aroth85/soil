import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import tasks


def create_optitype_workflow(bam_file, hla_type_file, is_rna=False, threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['optitype', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.commandline(
        name='extract_chr6',
        args=(
            'samtools', 'view', '-bh', '-f', '2', '-F', '4',
            mgd.InputFile(bam_file),
            '6',
            '|',
            'samtools', 'collate', '-O', '-', mgd.TempSpace('chr6_collate_temp'),
            '|',
            'samtools', 'bam2fq',
            '-1', mgd.TempOutputFile('chr6_reads_1.fq'),
            '-2', mgd.TempOutputFile('chr6_reads_2.fq'),
            '-',
        ),
    )

    workflow.transform(
        name='optitype',
        ctx={'mem': 24, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        func=tasks.run_optitype,
        args=(
            mgd.TempInputFile('chr6_reads_1.fq'),
            mgd.TempInputFile('chr6_reads_2.fq'),
            mgd.OutputFile(hla_type_file),
            mgd.TempSpace('optitype_temp'),
        ),
        kwargs={
            'is_rna': is_rna,
            'threads': threads,
        }
    )

    return workflow
