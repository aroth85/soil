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
            mgd.TempOutputFile('decoy.no_index.fasta'),
        ),
        kwargs={
            'decoy_prefix': 'XXX_'
        }
    )

    workflow.transform(
        name='index_decoy_db',
        func=tasks.build_index,
        args=(
            mgd.TempInputFile('decoy.no_index.fasta'),
            mgd.TempOutputFile('decoy.fasta')
        ),
        kwargs={
            'add_decoys': False
        }
    )

    workflow.transform(
        name='index_target_db',
        func=tasks.build_index,
        args=(
            mgd.InputFile(in_fasta_file),
            mgd.TempOutputFile('target.fasta')
        ),
        kwargs={
            'add_decoys': False
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
            mgd.TempInputFile('target.fasta'),
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
        name='convert_to_tsv_decoy',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.convert_mzid_to_tsv,
        args=(
            mgd.TempInputFile('decoy_search_results.mzid', 'split'),
            mgd.TempOutputFile('decoy_search.tsv', 'split'),
        )
    )

    workflow.transform(
        name='convert_to_tsv_target',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.convert_mzid_to_tsv,
        args=(
            mgd.TempInputFile('target_search_results.mzid', 'split'),
            mgd.TempOutputFile('target_search.tsv', 'split'),
        )
    )

    workflow.transform(
        name='merge_results',
        func=tasks.merge_results,
        args=(
            [mgd.TempInputFile('decoy_search.tsv', 'split'), mgd.TempInputFile('target_search.tsv', 'split')],
            mgd.TempOutputFile('merged.tsv')
        )
    )

    workflow.transform(
        name='convert_output',
        func=tasks.convert_msgf_to_final,
        args=(
            mgd.TempInputFile('merged.tsv'),
            mgd.OutputFile(out_file.replace('.tsv', '.msgf.tsv.gz'))
        )
    )

    workflow.transform(
        name='run_msgf2pin',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.percolator.tasks.convert_msgf_to_pin,
        args=(
            mgd.TempInputFile('decoy_search_results.mzid', 'split'),
            mgd.TempInputFile('target_search_results.mzid', 'split'),
            mgd.TempOutputFile('percolator_input.tsv'),
            mgd.TempSpace('msgf2pin_tmp')
        )
    )

    workflow.transform(
        name='run_percolator',
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.percolator.tasks.run_percolator,
        args=(
            mgd.TempInputFile('percolator_input.tsv'),
            mgd.TempOutputFile('final.tsv')
        )
    )

    workflow.transform(
        name='clean_up_decoy',
        func=tasks.clean_up,
        args=(
            [mgd.TempInputFile('decoy.fasta'), mgd.TempInputFile('target.fasta')],
            mgd.TempInputFile('final.tsv'),
            mgd.OutputFile(out_file)
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
        name='index_db',
        func=tasks.build_index,
        args=(
            mgd.InputFile(in_fasta_file),
            mgd.TempOutputFile('db.fasta')
        ),
        kwargs={
            'add_decoys': add_decoys
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

    workflow.transform(
        name='run_msgf_plus',
        axes=('split',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_search,
        args=(
            mgd.TempInputFile('db.fasta'),
            mgd.TempInputFile('spec_data.mzml', 'split'),
            mgd.TempOutputFile('search.mzid', 'split'),
            mgd.TempSpace('msgf_tmp', 'split'),
        ),
        kwargs={
            'add_decoys': False,
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
            mgd.TempOutputFile('merged.tsv')
        )
    )

    workflow.transform(
        name='convert_output',
        func=tasks.convert_msgf_to_final,
        args=(
            mgd.TempInputFile('merged.tsv'),
            mgd.TempOutputFile('final.tsv.gz')
        )
    )

    workflow.transform(
        name='clean_up',
        func=tasks.clean_up,
        args=(
            mgd.TempInputFile('db.fasta'),
            mgd.TempInputFile('final.tsv.gz'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow
