#!/usr/bin/env python
import pypeliner

import pot.workflows.wgs_align


def main(args):
    workflow = pot.workflows.wgs_align.create_single_lane_alignment_workflow(
        args.fastq_file_1,
        args.fastq_file_2,
        args.ref_genome_fasta_file,
        args.out_bam_file,
        bwa_threads=8,
        sambamba_threads=4
    )

    config = {
        'maxjobs': args.max_jobs,
        'nativespec': args.native_spec,
        'submit': args.submit,
        'tmpdir': args.working_dir,
    }

    pyp = pypeliner.app.Pypeline(config=config)

    pyp.run(workflow)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-f1', '--fastq_file_1', help='Path to 1st of paired end fastqs')

    parser.add_argument('-f2', '--fastq_file_2', help='Path to 1st of paired end fastqs')

    parser.add_argument('-r', '--ref_genome_fasta_file', help='Path to reference genome to align against')

    parser.add_argument('-o', '--out_bam_file', help='Path where output will be written in BAM format')

    parser.add_argument('--max_jobs', type=int, default=1)

    parser.add_argument('--native_spec', default='')

    parser.add_argument('--submit', choices=['drmaa', 'local'])

    parser.add_argument('--working_dir', default=None)

    cli_args = parser.parse_args()

    main(cli_args)
