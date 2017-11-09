'''
Created on 24 Apr 2017

@author: Andrew Roth
'''
import pypeliner
import pypeliner.managed as mgd
import pyteomics.fasta

import soil.utils.workflow
import soil.wrappers.percolator.tasks
import soil.wrappers.proteowizard.tasks

import tasks


def create_search_workflow(
        in_db_file,
        in_mzml_file,
        out_percolator_file,
        config={},
        split_size=1000):

    msgf_plus_config = config.get('msgf+', {})

    msgf_plus_config.update({
        'add_decoys': False,
        'add_features': True,
        'num_threads': 1,
    })

    sandbox = soil.utils.workflow.get_sandbox(['msgf_plus', 'percolator', 'proteowizard'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.subworkflow(
        name='build_db_index',
        func=create_index_workflow,
        args=(
            mgd.InputFile(in_db_file),
            mgd.TempOutputFile('decoy.fasta'),
            mgd.TempOutputFile('target.fasta'),
        )
    )

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

    workflow.transform(
        name='run_msgf_plus_decoy',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('decoy.fasta'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('decoy_search_results.mzid', 'split'),
            mgd.TempSpace('msgf_decoy_tmp', 'split'),
        ),
        kwargs=msgf_plus_config
    )

    workflow.transform(
        name='run_msgf_plus_target',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('target.fasta'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('target_search_results.mzid', 'split'),
            mgd.TempSpace('msgf_target_tmp', 'split'),
        ),
        kwargs=msgf_plus_config
    )

    workflow.transform(
        name='run_msgf2pin',
        axes=(),
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.percolator.tasks.convert_msgf_to_pin,
        args=(
            mgd.TempInputFile('decoy_search_results.mzid', 'split'),
            mgd.TempInputFile('target_search_results.mzid', 'split'),
            mgd.TempOutputFile('percolator_input.tsv'),
            mgd.TempSpace('msgf2pin_tmp'),
        )
    )

    workflow.commandline(
        name='concat_decoy_target_dbs',
        args=(
            'cat',
            mgd.TempInputFile('decoy.fasta'),
            mgd.TempInputFile('target.fasta'),
            '>',
            mgd.TempOutputFile('target_decoy.fasta'),
        )
    )

    workflow.transform(
        name='run_percolator',
        axes=(),
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.percolator.tasks.run_percolator,
        args=(
            mgd.TempInputFile('percolator_input.tsv'),
            mgd.OutputFile(out_percolator_file),
        ),
        kwargs={
            'db_file': mgd.TempInputFile('target_decoy.fasta'),
        }
    )

    return workflow


def create_index_workflow(db_file, decoy_db_file, target_db_file):
    """ Build separate target and decoy databases and index for MSGF+ search.
    """

    sandbox = soil.utils.workflow.get_sandbox(['msgf_plus'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.commandline(
        name='copy_db',
        args=(
            'cp',
            mgd.InputFile(db_file),
            mgd.OutputFile(target_db_file),
        )
    )

    workflow.transform(
        name='build_decoy_db',
        func=pyteomics.fasta.write_decoy_db,
        args=(
            mgd.InputFile(target_db_file),
            mgd.OutputFile(decoy_db_file),
        ),
        kwargs={
            'mode': 'reverse',
            'decoy_only': True,
            'file_mode': 'w',
        }
    )

    workflow.transform(
        name='index_decoy_db',
        func=tasks.build_index,
        args=(
            mgd.InputFile(decoy_db_file),
            mgd.OutputFile(decoy_db_file + '.index.done'),
        )
    )

    workflow.transform(
        name='index_target_db',
        func=tasks.build_index,
        args=(
            mgd.InputFile(target_db_file),
            mgd.OutputFile(target_db_file + '.index.done'),
        )
    )

    return workflow
