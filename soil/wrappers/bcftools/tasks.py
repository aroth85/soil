import os
import pypeliner.commandline as cli

from soil.wrappers.samtools.tasks import concatenate_vcf


def rename_chroms(chrom_map_file, in_file, out_file):
    if not os.path.exists(in_file + '.csi'):
        tmp_index = True

        cmd = ['bcftools', 'index', in_file]

        cli.execute(*cmd)

    else:
        tmp_index = False

    cmd = [
        'bcftools',
        'annotate',
        '--rename-chrs', chrom_map_file,
        in_file,
        '|',
        'bcftools',
        'view'
    ]

    out_file_type = _detect_file_type(out_file)

    if out_file_type == 'bcf':
        cmd.extend(['-O', 'b'])

    else:
        cmd.extend(['-O', 'z'])

    cmd.extend(['>', out_file])

    cli.execute(*cmd)

    if tmp_index:
        os.unlink(in_file + '.csi')


def _detect_file_type(file_name):
    if file_name.endswith('.bcf'):
        return 'bcf'

    elif file_name.endswith('.vcf.gz'):
        return 'vcf'

    else:
        return Exception('Unknown file extension: {}'.format(file_name))
