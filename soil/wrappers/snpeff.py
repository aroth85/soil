import os
import pypeliner.commandline as cli


def snpeff(db, in_vcf_file, out_file):
    os.environ['MALLOC_ARENA_MAX'] = '2'

    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-Xms2g',
        '-Xmx5g',
        '-hgvs1LetterAa',
        db,
        in_vcf_file,
        '>',
        out_file
    ]

    cli.execute(*cmd)
