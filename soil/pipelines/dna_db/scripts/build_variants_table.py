from collections import OrderedDict

import os
import pandas as pd
import shutil
import tempfile


def main(args):
    tmp_cache_dir = None

    try:
        if args.pyensembl_cache_dir is not None:
            tmp_cache_dir = tempfile.mkdtemp()

            shutil.rmtree(tmp_cache_dir)

            shutil.copytree(args.pyensembl_cache_dir, tmp_cache_dir)

            pyensembl_cache_dir = tmp_cache_dir

            os.environ['PYENSEMBL_CACHE_DIR'] = os.path.abspath(pyensembl_cache_dir)

        _build_variant_table(args.in_file, args.out_file, genome_version=args.genome_version)

    finally:
        if tmp_cache_dir is not None:
            shutil.rmtree(tmp_cache_dir)


def _build_variant_table(in_file, out_file, genome_version='GRCh37'):
    import varcode

    variants = varcode.load_vcf(in_file, genome=genome_version)

    effects = variants.effects()

    effects = effects.drop_silent_and_noncoding()

    df = []

    for eff in effects:
        if not eff.modifies_protein_sequence:
            continue

        row = OrderedDict((
            ('gene_id', eff.gene.gene_id),
            ('gene_name', eff.gene.gene_name),
            ('transcript_id', eff.transcript_id),
            ('transcript_name', eff.transcript_name),
            ('protein_id', eff.transcript.protein_id),
            ('chrom', 'chr{}'.format(eff.gene.contig)),
            ('nuc_variant', eff.variant.short_description),
            ('aa_variant', eff.short_description),
            ('beg', eff.variant.start),
            ('end', eff.variant.end),
            ('gene_beg', eff.gene.start),
            ('gene_end', eff.gene.end),
            ('nuc_ref', eff.variant.ref),
            ('nuc_alt', eff.variant.alt),
            ('aa_ref', getattr(eff, 'aa_ref', '')),
            ('aa_alt', getattr(eff, 'aa_alt', '')),
            ('aa_mutation_beg_offset', eff.aa_mutation_start_offset),
            ('aa_mutation_end_offset', eff.aa_mutation_end_offset),
            ('prot_ref', eff.original_protein_sequence),
            ('prot_alt', eff.mutant_protein_sequence),
            ('effect_type', str(eff).split("(")[0])
        ))

        df.append(row)

    df = pd.DataFrame(df)

    df['effect_type'] = df['effect_type'].str.lower()

    df = df.drop_duplicates(['protein_id', 'aa_variant'])

    df.to_csv(out_file, compression='gzip', index=False, sep='\t')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--in-file', required=True)

    parser.add_argument('-o', '--out-file', required=True)

    parser.add_argument('--genome-version', default='GRCh37')

    parser.add_argument('--pyensembl-cache-dir', default=None)

    args = parser.parse_args()

    main(args)
