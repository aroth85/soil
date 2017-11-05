#!/usr/bin/env python
import pypeliner

import soil.workflows.strelka_somatic


def main(args):
    workflow = soil.workflows.strelka_somatic.create_workflow(
        args.normal_bam_file,
        args.tumour_bam_file,
        args.ref_genome_fasta_file,
        args.out_vcf_file,
        chromosomes=args.chromosomes,
        split_size=int(1e7)
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

    parser.add_argument('-nb', '--normal_bam_file', required=True)

    parser.add_argument('-tb', '--tumour_bam_file', required=True)

    parser.add_argument('-r', '--ref_genome_fasta_file', required=True)

    parser.add_argument('-r', '--ref_genome_fasta_file', required=True)

    parser.add_argument('-o', '--out_vcf_file', required=True)

    parser.add_argument('--chromosomes', default=None, nargs='+')

    parser.add_argument('--max_jobs', type=int, default=1)

    parser.add_argument('--native_spec', default='')

    parser.add_argument('--submit', choices=['drmaa', 'local'])

    parser.add_argument('--working_dir', default=None)

    cli_args = parser.parse_args()

    main(cli_args)
