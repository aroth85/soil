import os
import pypeliner.commandline as cli


def call_variants(bam_file, ref_genome_fasta_file, log_file, out_file, region, rna=False):
    if not os.path.exists(bam_file + '.bai'):
        tmp_index = True

        cli.execute('samtools', 'index', bam_file)

    else:
        tmp_index = False

    cmd = [
        'platypus',
        'callVariants',
        '--bamFiles', bam_file,
        '--logFileName', log_file,
        '--refFile', ref_genome_fasta_file,
        '--regions', region,
        '-o', out_file,
    ]

    if rna:
        cmd.extend([
            '--filterDuplicates', 0,
            '--minMapQual', 0,
            '--minFlank', 0,
            '--maxReadLength', 500,
            '--minGoodQualBases', 10,
            '--minBaseQual', 20,
        ])

    cli.execute(*cmd)

    if tmp_index:
        os.unlink(bam_file + '.bai')
