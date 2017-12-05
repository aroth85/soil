'''
Created on 27 Jul 2017

@author: Andrew Roth
'''
import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow


def create_basic_workflow(fastq_file_1, fastq_file_2, out_file, threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['mixcr', ])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.commandline(
        name='align',
        ctx={'mem': 32, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        args=(
            'mixcr',
            'align',
            '-f',
            '-t', threads,
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            mgd.TempOutputFile('alignments.vdjca')
        )
    )

    workflow.commandline(
        name='assemble',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        args=(
            'mixcr',
            'assemble',
            '-f',
            '-t', 1,
            mgd.TempInputFile('alignments.vdjca'),
            mgd.TempOutputFile('clones.clns')
        )
    )

    workflow.commandline(
        name='export',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        args=(
            'mixcr',
            'exportClones',
            '-f',
            mgd.TempInputFile('clones.clns'),
            mgd.TempOutputFile('results.tsv')
        )
    )

    workflow.commandline(
        name='compress',
        args=(
            'gzip', '-c',
            mgd.TempInputFile('results.tsv'),
            '>',
            mgd.OutputFile(out_file)
        )
    )

    return workflow


def create_rnaseq_workflow(fastq_file_1, fastq_file_2, out_file, threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['mixcr', ])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.commandline(
        name='align',
        ctx={'mem': 32, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        args=(
            'mixcr',
            'align',
            '-p', 'rna-seq',
            '-s', 'hsa',
            '-OallowPartialAlignments=true',
            '-f',
            '-t', threads,
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            mgd.TempOutputFile('alignments.vdjca')
        )
    )

    workflow.commandline(
        name='assemblePartial_1',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        args=(
            'mixcr',
            'assemblePartial',
            '-f',
            mgd.TempInputFile('alignments.vdjca'),
            mgd.TempOutputFile('alignments_rescued_1.vdjca')
        )
    )

    workflow.commandline(
        name='assemblePartial_2',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        args=(
            'mixcr',
            'assemblePartial',
            '-f',
            mgd.TempInputFile('alignments_rescued_1.vdjca'),
            mgd.TempOutputFile('alignments_rescued_2.vdjca')
        )
    )

    workflow.commandline(
        name='extendAlignments',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        args=(
            'mixcr',
            'extendAlignments',
            '-f',
            mgd.TempInputFile('alignments_rescued_2.vdjca'),
            mgd.TempOutputFile('alignments_rescued_2_extended.vdjca')
        )
    )

    workflow.commandline(
        name='assemble',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        args=(
            'mixcr',
            'assemble',
            '-f',
            '-t', threads,
            mgd.TempInputFile('alignments_rescued_2_extended.vdjca'),
            mgd.TempOutputFile('clones.clns')
        )
    )

    workflow.commandline(
        name='export',
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        args=(
            'mixcr',
            'exportClones',
            '-f',
            mgd.TempInputFile('clones.clns'),
            mgd.TempOutputFile('results.tsv')
        )
    )

    workflow.commandline(
        name='compress',
        args=(
            'gzip', '-c',
            mgd.TempInputFile('results.tsv'),
            '>',
            mgd.OutputFile(out_file)
        )
    )

    return workflow
