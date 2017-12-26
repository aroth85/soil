import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow

import tasks


def create_topiary_workflow(hla_alleles, in_file, out_file, iedb_dir=None, genome='GRCh37', pyensembl_cache_dir=None):
    """ Run topiary.

    Parameters
    ----------
    hla_alleles: list
        List of HLA alleles i.e. A*02:01.
    in_file: str
        Path to VCF file with variants.
    out_file: str
        Path where output will be written in tsv format.
    """
    sandbox = soil.utils.workflow.get_sandbox(['topiary', ])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(obj=mgd.TempOutputObj('hla_alleles'), value=hla_alleles)

    workflow.transform(
        name='run_topiary',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.run_topiary,
        args=(
            mgd.TempInputObj('hla_alleles'),
            mgd.InputFile(in_file),
            mgd.TempOutputFile('raw.tsv')
        ),
        kwargs={
            'iedb_dir': iedb_dir,
            'genome': genome,
            'predictor': 'netmhc',
            'pyensembl_cache_dir': pyensembl_cache_dir
        }
    )

    workflow.transform(
        name='reformat_output',
        func=tasks.reformat_output,
        args=(
            mgd.TempInputFile('raw.tsv'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow
