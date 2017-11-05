#!/usr/bin/env python
import soil.utils.cli
import soil.workflows.wgs_align


def main(args):
    workflow = soil.workflows.wgs_align.create_single_lane_alignment_workflow(
        args.fastq_file_1,
        args.fastq_file_2,
        args.ref_genome_fasta_file,
        args.out_bam_file,
        bwa_threads=8,
        sambamba_threads=4
    )

    soil.utils.cli.run_workflow(args, workflow)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-f1', '--fastq_file_1', help='Path to 1st of paired end fastqs', required=True)

    parser.add_argument('-f2', '--fastq_file_2', help='Path to 1st of paired end fastqs', required=True)

    parser.add_argument(
        '-r', '--ref_genome_fasta_file', help='Path to reference genome to align against', required=True)

    parser.add_argument('-o', '--out_bam_file', help='Path where output will be written in BAM format', required=True)

    soil.utils.cli.add_pipeline_args(parser)

    cli_args = parser.parse_args()

    main(cli_args)
