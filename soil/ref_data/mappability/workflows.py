import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow

import tasks


def create_mappability_workflow(
        ref_genome_fasta_file,
        out_file,
        k=100,
        max_map_qual=None,
        split_size=int(1e7),
        threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['bwa', 'samtools', 'ucsc-bedgraphtobigwig'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='split_fasta_by_chrom',
        func=tasks.split_fasta_by_chrom,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('chrom.fasta', 'chrom'),
        )
    )

    workflow.transform(
        name='create_kmer_reads',
        axes=('chrom',),
        ctx={'mem': 4, 'mem_retry_increment': 2, 'num_retry': 3},
        func=tasks.create_kmer_reads,
        args=(
            mgd.TempInputFile('chrom.fasta', 'chrom'),
            mgd.TempOutputFile('reads.fa', 'chrom', 'kmer_group'),
        ),
        kwargs={
            'k': k,
            'split_size': split_size,
        }
    )

    workflow.transform(
        name='align_kmers',
        axes=('chrom', 'kmer_group'),
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        func=tasks.bwa_mem_align,
        args=(
            mgd.TempInputFile('reads.fa', 'chrom', 'kmer_group'),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('aligned.bam', 'chrom', 'kmer_group'),
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.transform(
        name='compute_mappability',
        axes=('chrom', 'kmer_group'),
        ctx={'mem': 4, 'mem_retry_increment': 2, 'num_retry': 3},
        func=tasks.compute_mappability,
        args=(
            mgd.TempInputFile('aligned.bam', 'chrom', 'kmer_group'),
            mgd.TempOutputFile('mappability.tsv', 'chrom', 'kmer_group'),
        ),
        kwargs={
            'max_map_qual': max_map_qual,
        }
    )

    workflow.transform(
        name='compute_chrom_mean_mappability',
        axes=('chrom',),
        ctx={'mem': 16, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.compute_chrom_mean_mappability,
        args=(
            mgd.TempInputFile('mappability.tsv', 'chrom', 'kmer_group'),
            mgd.TempOutputFile('mean_mappability.tsv', 'chrom'),
        ),
    )

    workflow.transform(
        name='write_bed',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.write_bed,
        args=(
            mgd.TempInputFile('mean_mappability.tsv', 'chrom'),
            mgd.TempOutputFile('mean_mappability.bed'),
        ),
    )

    workflow.transform(
        name='write_chrom_sizes',
        func=tasks.write_chrom_sizes,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('chrom_sizes.txt'),
        )
    )

    workflow.commandline(
        name='write_big_wig',
        args=(
            'bedGraphToBigWig',
            mgd.TempInputFile('mean_mappability.bed'),
            mgd.TempInputFile('chrom_sizes.txt'),
            mgd.OutputFile(out_file),
        )
    )

    return workflow
