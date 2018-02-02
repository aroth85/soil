'''
Created on 24 Apr 2017

@author: Andrew Roth
'''
import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.proteowizard.tasks

import tasks


def create_search_workflow(
        in_fasta_file,
        in_mzml_file,
        out_file,
        add_decoys=True,
        fixed_mods=None,
        max_mods=1,
        split_size=1000,
        variable_mods=None):

    sandbox = soil.utils.workflow.get_sandbox(['msgf_plus', 'proteowizard'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='split_mzml_file',
        axes=(),
        func=soil.wrappers.proteowizard.tasks.split_mzml_file,
        args=(
            mgd.InputFile(in_mzml_file),
            mgd.TempOutputFile('spec_data.mzml', 'split'),
            mgd.TempSpace('split_tmp'),
        ),
        kwargs={
            'split_size': split_size,
        }
    )

    workflow.commandline(
        name='copy_db',
        ctx={'local': True},
        axes=('split',),
        args=(
            'cp',
            mgd.InputFile(in_fasta_file),
            mgd.TempOutputFile('db.fasta', 'split'),
        )
    )

    workflow.transform(
        name='run_msgf_plus',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('db.fasta', 'split'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('search.mzid', 'split'),
            mgd.TempSpace('msgf_tmp', 'split'),
        ),
        kwargs={
            'add_decoys': add_decoys,
            'fixed_mods': fixed_mods,
            'max_mods': max_mods,
            'variable_mods': variable_mods
        }
    )

    workflow.transform(
        name='convert_to_tsv',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.convert_mzid_to_tsv,
        args=(
            mgd.TempInputFile('search.mzid', 'split'),
            mgd.TempOutputFile('search.tsv', 'split'),
        )
    )

    workflow.transform(
        name='merge_results',
        func=tasks.merge_results,
        args=(
            mgd.TempInputFile('search.tsv', 'split'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow
