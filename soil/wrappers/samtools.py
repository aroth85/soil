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
    in_files = soil.utils.workflow.flatten_input(in_files)

    tmp_index_files = []

    if allow_overlap:
        cmd = ['bcftools', 'concat', '-a', '-O', 'z', '-o', out_file]

        # bcftools requires that the input files be indexed when overlap is allowed
        for file_name in in_files:
            tmp_index_file_name = file_name + '.tbi'

            if not os.path.exists(tmp_index_file_name):
                index_vcf(file_name, index_file=tmp_index_file_name)

                tmp_index_files.append(tmp_index_file_name)

    else:
        cmd = ['bcftools', 'concat', '-O', 'z', '-o', out_file]

    cmd.extend(in_files)

    cli.execute(*cmd)

    for file_name in tmp_index_files:
        os.unlink(tmp_index_files)

    if index_file is not None:
        index_vcf(out_file, index_file=index_file)


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
