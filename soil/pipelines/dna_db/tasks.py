import pandas as pd
import pypeliner.commandline as cli

import soil.utils.package_data


def build_variant_fasta(in_file, out_file):
    df = pd.read_csv(in_file, sep='\t')

    df['prot_alt'] = df['prot_alt'].fillna('')

    df['prot_alt'] = df['prot_alt'].astype(str)

    with open(out_file, 'w') as out_fh:
        for _, row in df.iterrows():
            if len(row['prot_alt']) == 0:
                continue

            header = 'mut|{0}_{1}'.format(row['protein_id'], row['aa_variant'])

            out_fh.write('>{}\n'.format(header))

            out_fh.write(row['prot_alt'] + '\n')


def build_variant_table(in_file, out_file, genome_version='GRCh37', pyensembl_cache_dir=None):
    script = soil.utils.package_data.load_data_file('pipelines/dna_db/scripts/build_variants_table.py')

    cmd = [
        'python',
        script,
        '-i', in_file,
        '-o', out_file,
        '--genome-version', genome_version
    ]

    if pyensembl_cache_dir is not None:
        cmd.extend(['--pyensembl-cache-dir', pyensembl_cache_dir])

    cli.execute(*cmd)


def clean_ref_proteome_ids(in_file, out_file):
    with open(in_file) as in_fh, open(out_file, 'w') as out_fh:
        for line in in_fh:
            line = line.strip()

            if line.startswith('>'):
                line = '>' + '|'.join(['ref', line[1:].split('|')[0].split('.')[0]])

            out_fh.write(line + '\n')
