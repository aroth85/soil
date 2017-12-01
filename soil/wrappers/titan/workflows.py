import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.platypus.workflows

import tasks


def create_titan_workflow(
        normal_bam_file,
        tumour_bam_file,
        dbsnp_vcf_file,
        mappability_file,
        ref_genome_fasta_file,
        out_sentinel_file,
        out_stats_file,
        exome_bed_file=None,
        sample='Tumour',
        out_normal_vcf_file=None):

    if out_normal_vcf_file is None:
        normal_vcf = mgd.TempFile('normal.vcf.gz')

    else:
        normal_vcf = mgd.File(out_normal_vcf_file)

    sandbox = soil.utils.workflow.get_sandbox(['hmmcopy', 'hmmcopy_utils', 'snpsift', 'titan'])

    sandbox.packages.extend(['pandas', 'rpy2'])

    chromosomes = soil.utils.genome.get_autosomal_chromosomes(normal_bam_file)

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(
        obj=mgd.TempOutputObj('init_params', 'param_idx'),
        value=tasks.create_intialization_parameters()
    )

    workflow.subworkflow(
        name='call_snps',
        func=soil.wrappers.platypus.workflows.create_single_sample_workflow,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            normal_vcf.as_output(),
        ),
        kwargs={
            'chromosomes': chromosomes,
            'split_size': int(1e7),
        }
    )

    workflow.commandline(
        name='annotate_dbsnp_status',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'annotate',
            mgd.InputFile(dbsnp_vcf_file),
            normal_vcf.as_input(),
            '>',
            mgd.TempOutputFile('normal.dbsnp.vcf'),
        ),
    )

    workflow.commandline(
        name='annotate_variant_type',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'varType',
            mgd.TempInputFile('normal.dbsnp.vcf'),
            '>',
            mgd.TempOutputFile('normal.dbsnp.vartype.vcf'),
        ),
    )

    workflow.commandline(
        name='filter_het_snps',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'filter',
            "isHet(GEN[0]) & ((exists ID) & ( ID =~ 'rs' )) & (exists SNP)",
            mgd.TempInputFile('normal.dbsnp.vartype.vcf'),
            '>',
            mgd.TempOutputFile('het.snps.vcf'),
        ),
    )

    workflow.transform(
        name='split_vcf',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.split_vcf,
        args=(
            mgd.TempInputFile('het.snps.vcf'),
            mgd.TempOutputFile('split.vcf', 'split'),
            mgd.TempSpace('split_tmp')
        ),
        kwargs={
            'split_size': int(1e4),
        }
    )

    workflow.transform(
        name='get_allele_counts',
        axes=('split',),
        func=tasks.get_snv_allele_counts_for_vcf_targets,
        args=(
            mgd.InputFile(tumour_bam_file),
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('split.tsv', 'split'),
        )
    )

    workflow.transform(
        name='merge_counts',
        func=tasks.merge_counts,
        args=(
            mgd.TempInputFile('split.tsv', 'split'),
            mgd.TempOutputFile('allele_counts.tsv'),
        )
    )

    workflow.commandline(
        name='build_normal_wig',
        args=(
            'readCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(normal_bam_file),
            '>',
            mgd.TempOutputFile('normal.wig'),
        )
    )

    workflow.commandline(
        name='build_tumour_wig',
        args=(
            'readCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(tumour_bam_file),
            '>',
            mgd.TempOutputFile('tumour.wig'),
        )
    )

    workflow.commandline(
        name='build_gc_wig',
        args=(
            'gcCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(ref_genome_fasta_file),
            '>',
            mgd.TempOutputFile('gc.wig'),
        )
    )

    workflow.commandline(
        name='build_mappability_wig',
        args=(
            'mapCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(mappability_file),
            '>',
            mgd.TempOutputFile('mappability.wig'),
        )
    )

    workflow.transform(
        name='build_coverage_file',
        func=tasks.build_coverage_file,
        args=(
            mgd.TempInputFile('normal.wig'),
            mgd.TempInputFile('tumour.wig'),
            mgd.TempInputFile('gc.wig'),
            mgd.TempInputFile('mappability.wig'),
            mgd.TempOutputFile('coverage.wig'),
        ),
        kwargs={
            'target_file': exome_bed_file,
        }
    )

    workflow.transform(
        name='run_titan',
        axes=('param_idx',),
        func=tasks.run_titan,
        args=(
            mgd.TempInputFile('coverage.wig'),
            mgd.TempInputFile('allele_counts.tsv'),
            mgd.TempInputObj('init_params', 'param_idx'),
            mgd.TempOutputFile('run.tar.gz', 'param_idx'),
            mgd.TempSpace('titan_tmp', 'param_idx'),
        ),
        kwargs={
            'is_exome': (exome_bed_file is not None),
            'sample': sample,
        }
    )

    workflow.transform(
        name='build_run_stats_file',
        func=tasks.build_run_stats_file,
        args=(
            mgd.TempInputFile('run.tar.gz', 'param_idx'),
            mgd.TempInputObj('init_params', 'param_idx'),
            mgd.OutputFile(out_stats_file),
        ),
    )

    workflow.transform(
        name='select_solution',
        func=tasks.select_optimal_run,
        args=(
            mgd.InputFile(out_stats_file),
        ),
        ret=mgd.TempOutputObj('optimal_idx'),
    )

    workflow.transform(
        name='copy_optimal_solution',
        func=tasks.copy_optimal_solution,
        args=(
            mgd.TempInputObj('optimal_idx'),
            mgd.TempInputFile('run.tar.gz', 'param_idx'),
            mgd.OutputFile(out_sentinel_file),
            #             mgd.OutputFile(out_params_file),
            #             mgd.OutputFile(out_segs_file),
            #             mgd.OutputFile(out_plots_dir),
        )
    )

#     workflow.transform(
#         name='run_titan_with_optimal_params',
#         func=tasks.run_titan,
#         args=(
#             mgd.TempInputFile('coverage.wig'),
#             mgd.TempInputFile('allele_counts.tsv'),
#             mgd.TempInputObj('optimal_params'),
#             mgd.OutputFile(out_params_file),
#             mgd.OutputFile(out_segs_file),
#             mgd.TempSpace('titan_tmp'),
#         ),
#         kwargs={
#             'estimate_normal_contam': True,
#             'estimate_ploidy': True,
#             'is_exome': (exome_bed_file is not None),
#             'plot_dir': mgd.OutputFile(out_plots_dir),
#             'sample': sample,
#         }
#     )

    return workflow