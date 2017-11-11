import Bio.Seq
import Bio.SeqIO
import pypeliner.commandline as cli
import pysftp


def decompress(in_file, out_file):
    if in_file.endswith('gz'):
        cmd = ['gzip', '-cd', in_file, '>', out_file]

    cli.execute(*cmd)


def download(url, local_path):
    cli.execute('wget', url, '-c', '-O', local_path)


def download_from_sftp(host, host_path, local_path, user, password):
    with pysftp.Connection(host, username=user, password=password) as sftp:
        sftp.get(host_path, localpath=local_path)


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
