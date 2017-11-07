import os
import pypeliner.commandline as cli


def mpileup2snp(in_file, out_file):
    if os.stat(in_file).st_size == 0:
        cmd = ['echo', 'a']

    else:
        cmd = ['cat', in_file]

    cmd.extend([
        '|'
        'varscan',
        'mpileup2snp',
        '--output-vcf', 1,
        '>',
        out_file
    ])

    cli.execute(*cmd)
