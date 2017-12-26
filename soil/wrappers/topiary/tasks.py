import os
import pandas as pd
import pypeliner.commandline as cli
import soil.utils.file_system


def run_topiary(
        hla_alleles,
        in_file,
        out_file,
        iedb_dir=None,
        genome='GRCh37',
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
        '--mhc-predictor', predictor
    ]

    cli.execute(*cmd)


def reformat_output(in_file, out_file):
    df = pd.read_csv(in_file, sep=',')

    df = df.drop('#', axis=1)

    df.to_csv(out_file, index=False, sep='\t')
