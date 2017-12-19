import pysam
import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow

import tasks

default_ctx = {'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3}


def create_ref_panel_phase_workflow(genetic_map_file, ref_file, target_file, out_file):
    """ Run EAGLE using a reference panel.
    """

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'eagle-phase'])

    workflow = pypeliner.workflow.Workflow(default_ctx=default_ctx, default_sandbox=sandbox)

    workflow.setobj(
        obj=mgd.TempOutputObj('chrom', 'chrom'),
        value=get_chromosomes(target_file)
    )

    workflow.transform(
        name='split_ref',
        axes=('chrom',),
        func=tasks.get_chrom_variant_file,
        args=(
            mgd.TempInputObj('chrom', 'chrom'),
            mgd.InputFile(ref_file),
            mgd.TempOutputFile('ref.bcf', 'chrom')
        )
    )

    workflow.transform(
        name='split_target',
        axes=('chrom',),
        func=tasks.get_chrom_variant_file,
        args=(
            mgd.TempInputObj('chrom', 'chrom'),
            mgd.InputFile(target_file),
            mgd.TempOutputFile('target.bcf', 'chrom')
        )
    )

    workflow.transform(
        name='run_eagle',
        axes=('chrom',),
        func=tasks.run_eagle,
        args=(
            mgd.InputFile(genetic_map_file),
            mgd.TempInputFile('ref.bcf', 'chrom'),
            mgd.TempInputFile('target.bcf', 'chrom'),
            mgd.TempOutputFile('phased.bcf', 'chrom'),
            mgd.TempSpace('eagle_tmp', 'chrom')
        )
    )

    workflow.transform(
        name='concat_results',
        func=tasks.concat_results,
        args=(
            mgd.TempInputFile('phased.bcf', 'chrom'),
            mgd.OutputFile(out_file)
        )
    )

    workflow.commandline(
        name='index',
        args=(
            'bcftools',
            'index',
            '-o', mgd.OutputFile(out_file + '.csi'),
             mgd.InputFile(out_file)
        )
    )

    return workflow


def get_chromosomes(variant_file):
    vf = pysam.VariantFile(variant_file, 'r')

    chroms = []

    autosomes = [str(x) for x in range(1, 23)]

    for raw_chrom in list(vf.header.contigs):
        parsed_chrom = raw_chrom.replace('chr', '')

        if parsed_chrom in autosomes:
            chroms.append(raw_chrom)

    return chroms
