"""
Wrappers for the samtools suite of tools which include bcftools, smatools and tabix.

http://www.htslib.org
"""
import os
import pypeliner.commandline as cli
import shutil

import soil.utils.workflow


def concatenate_vcf(in_files, out_file, allow_overlap=False, index_file=None):
    """ Concatenate VCF files.

    :param in_files: dict with values being files to be concatenated. Files will be concatenated based on sorted order
        of keys.
    :param out_file: path where output file will be written in VCF format.
    :param allow_overlap: bool indicating whether there maybe overlapping regions in the files.
    :param index_file: path where tabix file will be written.
    """
    if allow_overlap:
        cmd = ['bcftools', 'concat', '-a', '-O', 'z', '-o', out_file]

    else:
        cmd = ['bcftools', 'concat', '-O', 'z', '-o', out_file]

    tmp_index_files = []

    for file_name in in_files:
        if not os.path.exists(file_name + '.tbi'):
            tmp_index_files.append(file_name + '.tbi')

            cli.execute('tabix', '-f', '-p', 'vcf', file_name)

    cmd += soil.utils.workflow.flatten_input(in_files)

    cli.execute(*cmd)

    if index_file is not None:
        index_vcf(out_file, index_file=index_file)

    for file_name in tmp_index_files:
        os.unlink(file_name)


def index_vcf(vcf_file, index_file=None):
    """ Create a tabix index for a VCF file

    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.

    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.

    """
    cli.execute('tabix', '-f', '-p', 'vcf', vcf_file)

    if index_file is None:
        _rename_index(vcf_file, '.tbi')

    else:
        shutil.move(vcf_file + '.tbi', index_file)


def _rename_index(in_file, index_suffix):
    if in_file.endswith('.tmp'):
        index_file = in_file[:-4] + index_suffix

        try:
            os.remove(index_file)
        except:
            pass

        os.rename(in_file + index_suffix, index_file)
