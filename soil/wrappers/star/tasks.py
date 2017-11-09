import os
import pypeliner.commandline as cli
import shutil


def align(
        fastq_file_1,
        fastq_file_2,
        ref_genome_dir,
        out_file,
        tmp_dir,
        log_dir=None,
        read_group_info=None,
        threads=1,
        unaligned_read_fastq_1=None,
        unaligned_read_fastq_2=None):
    """ Align paired end reads with the STAR aligner
    """

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_prefix = os.path.join(tmp_dir, 'alignment')

    cmd = [
        'STAR',
        '--runThreadN', threads,
        '--genomeDir', ref_genome_dir,
        '--readFilesIn', fastq_file_1, fastq_file_2,
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outSAMattributes', 'NH', 'HI', 'NM', 'MD', 'AS', 'nM',
        '--readFilesCommand', 'zcat',
        '--outFileNamePrefix', tmp_prefix,
    ]

    if read_group_info is not None:
        read_group_str = ['ID:{0}'.format(read_group_info['ID']), ]
        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue

            read_group_str.append(':'.join((key, value)))

        read_group_str = '\t'.join(read_group_str)

        cmd.extend(['--outSAMattrRGline', read_group_str])

    if unaligned_read_fastq_1 is not None:
        cmd.extend(['--outReadsUnmapped', 'Fastx'])

    cli.execute(*cmd)

    tmp_out_file = tmp_prefix + 'Aligned.sortedByCoord.out.bam'

    shutil.move(tmp_out_file, out_file)

    if unaligned_read_fastq_1 is not None:
        tmp_out_file = tmp_prefix + 'Unmapped.out.mate1'

        shutil.move(tmp_out_file, unaligned_read_fastq_1)

    if unaligned_read_fastq_2 is not None:
        tmp_out_file = tmp_prefix + 'Unmapped.out.mate2'

        shutil.move(tmp_out_file, unaligned_read_fastq_2)

    if log_dir is not None:
        if os.path.exists(log_dir):
            shutil.rmtree(log_dir)

        shutil.copytree(tmp_dir, log_dir)


def index(ref_genome_fasta_file, transcript_gtf_file, out_sentinel_file, overhang=100, threads=1):
    """ Build an index of a FASTA file for use with the STAR alignment programs.
    """

    cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--runThreadN', threads,
        '--genomeDir', os.path.dirname(out_sentinel_file),
        '--genomeFastaFiles', ref_genome_fasta_file,
        '--sjdbGTFfile', transcript_gtf_file,
        '--sjdbOverhang', overhang,
    ]

    cli.execute(*cmd)

    open(out_sentinel_file, 'w').close()
