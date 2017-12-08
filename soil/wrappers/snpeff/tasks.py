import os
import pypeliner.commandline as cli


def snpeff(db, in_vcf_file, out_file, data_dir=None):
    os.environ['MALLOC_ARENA_MAX'] = '2'

    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-Xms2g',
        '-Xmx5g',
        '-hgvs1LetterAa',
    ]

    if data_dir is not None:
        cmd.extend(['-dataDir', data_dir])
        cmd.append('-nodownload')

    cmd.extend([
        db,
        in_vcf_file,
        '>',
        out_file
    ])

    cli.execute(*cmd)
