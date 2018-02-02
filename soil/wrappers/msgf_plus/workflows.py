import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.percolator.tasks
import soil.wrappers.proteowizard.tasks

import tasks


def create_percolator_workflow(
        in_fasta_file,
        in_mzml_file,
        out_file,
        fixed_mods=None,
        max_mods=1,
        split_size=1000,
        variable_mods=None):

    sandbox = soil.utils.workflow.get_sandbox(['msgf_plus', 'percolator', 'proteowizard'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='build_decoy_db',
        func=tasks.build_decoy_db,
        args=(
            mgd.InputFile(in_fasta_file),
            mgd.TempOutputFile('decoy.fasta'),
        ),
        kwargs={
            'decoy_prefix': 'XXX_'
        }
    )

    workflow.transform(
        name='split_mzml_file',
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
        name='copy_decoy_db',
        ctx={'local': True},
        axes=('split',),
        args=(
            'cp',
            mgd.TempInputFile('decoy.fasta'),
            mgd.TempOutputFile('decoy.fasta', 'split'),
        )
    )

    workflow.commandline(
        name='copy_target_db',
        ctx={'local': True},
        axes=('split',),
        args=(
            'cp',
            mgd.InputFile(in_fasta_file),
            mgd.TempOutputFile('target.fasta', 'split'),
        )
    )

    workflow.transform(
        name='run_msgf_plus_decoy',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('decoy.fasta', 'split'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('decoy_search_results.mzid', 'split'),
            mgd.TempSpace('msgf_decoy_tmp', 'split'),
        ),
        kwargs={
            'add_decoys': False,
            'add_features': True,
            'fixed_mods': fixed_mods,
            'max_mods': max_mods,
            'variable_mods': variable_mods
        }
    )

    workflow.transform(
        name='run_msgf_plus_target',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('target.fasta', 'split'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('target_search_results.mzid', 'split'),
            mgd.TempSpace('msgf_target_tmp', 'split'),
        ),
        kwargs={
            'add_decoys': False,
            'add_features': True,
            'fixed_mods': fixed_mods,
            'max_mods': max_mods,
            'variable_mods': variable_mods
        }
    )

    workflow.transform(
        name='run_msgf2pin',
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
            mgd.InputFile(in_fasta_file),
            mgd.TempInputFile('decoy.fasta'),
            '>',
            mgd.TempOutputFile('target_decoy.fasta'),
        )
    )

    workflow.transform(
        name='run_percolator',
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.percolator.tasks.run_percolator,
        args=(
            mgd.TempInputFile('percolator_input.tsv'),
            mgd.OutputFile(out_file),
        )
    )

    return workflow


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
