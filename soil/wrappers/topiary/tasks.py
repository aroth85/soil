import os
import pandas as pd
import pypeliner.commandline as cli
import soil.utils.file_system
import soil.utils.workflow


def run_topiary(
        hla_alleles,
        in_file,
        out_file,
        iedb_dir=None,
        genome='GRCh37',
        peptide_length=9,
        predictor='nethmhc',
        pyensembl_cache_dir=None):

    if iedb_dir is not None:
        if predictor == 'netmhc':
            exe_dir = os.path.dirname(soil.utils.file_system.find('netMHC-4.0.readme', iedb_dir))

        os.environ['PATH'] = ':'.join([exe_dir, os.environ['PATH']])

    if pyensembl_cache_dir is not None:
        os.environ['PYENSEMBL_CACHE_DIR'] = os.path.abspath(pyensembl_cache_dir)

    hla_alleles = ','.join(['HLA-{}'.format(x) for x in hla_alleles])

    cmd = [
        'topiary',

        '--mhc-alleles', hla_alleles,
        '--vcf', in_file,
        '--output-csv', out_file,

        '--genome', genome,
        '--mhc-predictor', predictor,
        '--mhc-peptide-lengths', peptide_length
    ]

    cli.execute(*cmd)


def reformat_output(in_files, out_file):
    data = []

    for file_name in soil.utils.workflow.flatten_input(in_files):
        df = pd.read_csv(file_name, sep=',')

        df = df.drop('#', axis=1)

        data.append(df)

    data = pd.concat(data)

    data.to_csv(out_file, index=False, sep='\t')
