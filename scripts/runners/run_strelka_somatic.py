#!/usr/bin/env python
import soil.utils.cli
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

    soil.utils.cli.run_workflow(args, workflow)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-nb', '--normal_bam_file', required=True)

    parser.add_argument('-tb', '--tumour_bam_file', required=True)

    parser.add_argument('-r', '--ref_genome_fasta_file', required=True)

    parser.add_argument('-r', '--ref_genome_fasta_file', required=True)

    parser.add_argument('-o', '--out_vcf_file', required=True)

    parser.add_argument('--chromosomes', default=None, nargs='+')

    soil.utils.cli.add_pipeline_args(parser)

    cli_args = parser.parse_args()

    main(cli_args)
