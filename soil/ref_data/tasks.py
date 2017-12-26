import Bio.Seq
import Bio.SeqIO
import os
import pypeliner.commandline as cli
import pysam
import pysftp
import shutil

import soil.ref_data.mappability.workflows


def configure_iedb_module(in_sentinel, out_sentinel):
    iedb_dir = os.path.dirname(in_sentinel)

    os.chdir(iedb_dir)

    cli.execute('tcsh', 'configure')

    open(out_sentinel, 'w').close()


def decompress(in_file, out_file):
    file_type = _guess_file_type(in_file)

    if file_type == 'gz':
        cmd = ['gzip', '-cd', in_file, '>', out_file]

    else:
        cmd = ['cp', in_file, out_file]

    cli.execute(*cmd)


def download(url, local_path):
    cli.execute('wget', url, '-c', '-O', local_path)

    cli.execute('touch', local_path)


def download_from_sftp(host, host_path, local_path, user, password):
    with pysftp.Connection(host, username=user, password=password) as sftp:
        sftp.get(host_path, localpath=local_path)


def download_pyensembl_cache(ensembl_version, out_sentinel):
    os.environ['PYENSEMBL_CACHE_DIR'] = os.path.dirname(out_sentinel)

    cmd = [
        'pyensembl',
        'install',
        '--release', ensembl_version,
        '--species', 'human'
    ]

    cli.execute(*cmd)

    open(out_sentinel, 'w').close()


def extract_tar_file(in_file, out_sentinel):
    out_dir = os.path.dirname(out_sentinel)

    cmd = [
        'tar',
        '-C', out_dir,
        '-zxvf', in_file,
        '--strip-components=1'
    ]

    cli.execute(*cmd)

    open(out_sentinel, 'w').close()


def filter_bad_proiteins(in_file, out_file):
    bad_chars = set(Bio.Alphabet.IUPAC.ExtendedIUPACProtein.letters) - set(Bio.Alphabet.IUPAC.IUPACProtein.letters)

    with open(in_file, 'r') as in_fh, open(out_file, 'w') as out_fh:
        for record in Bio.SeqIO.parse(in_fh, 'fasta'):
            write_record = True

            for aa in bad_chars:
                if aa in record.seq:
                    write_record = False

            if write_record:
                Bio.SeqIO.write(record, out_fh, 'fasta')


def lex_sort_fasta(in_file, out_file):
    raw_ref = pysam.FastaFile(in_file)

    lex_order = ['chr{}'.format(i) for i in range(1, 23) + ['X', 'Y', 'M']]

    contigs = set(raw_ref.references) - set(lex_order)

    lex_order = lex_order + sorted(contigs)

    with open(out_file, 'w') as out_fh:
        for i, chrom in enumerate(lex_order):
            record = Bio.SeqIO.SeqRecord(Bio.Seq.Seq(raw_ref.fetch(chrom)), description='', id=chrom)

            Bio.SeqIO.write(record, out_fh, 'fasta')

    lex_ref = pysam.FastaFile(out_file)

    for chrom in raw_ref.references:
        assert raw_ref.fetch(chrom) == lex_ref.fetch(chrom)


def mappability_wrapper(bwa_sentinel_file, out_file, **kwargs):
    """ Simple wrapper so task can depend on sentinel and ensure index is built
    """
    return soil.ref_data.mappability.workflows.create_mappability_workflow(
        bwa_sentinel_file.replace('.bwa_index.done', ''), out_file, **kwargs
    )


def unzip_file(in_file, out_sentinel, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    cmd = ['unzip', in_file, '-d', tmp_dir]

    cli.execute(*cmd)

    out_dir = os.path.dirname(out_sentinel)

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    shutil.copytree(tmp_dir, out_dir)

    open(out_sentinel, 'w').close()

    shutil.rmtree(tmp_dir)


def _guess_file_type(filename):
    magic_dict = {
        '\x1f\x8b\x08': 'gz',
        '\x42\x5a\x68': 'bz2',
        '\x50\x4b\x03\x04': 'zip'
    }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        file_start = f.read(max_len)

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return None
