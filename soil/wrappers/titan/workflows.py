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
        out_file,
        exome_bed_file=None,
        sample='Tumour',
        threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['hmmcopy', 'hmmcopy_utils', 'titan'])

    sandbox.channels.append(['conda-forge'])

    sandbox.packages.extend(['pandas', 'rpy2'])

    chromosomes = soil.utils.genome.load_bam_chromosome_lengths(normal_bam_file, 'autosomes')

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(
        obj=mgd.TempOutputObj('init_params', 'param_idx'),
        value=tasks.create_intialization_parameters()
    )

    workflow.subworkflow(
        name='get_allele_counts',
        func=create_allele_counts_workflow,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(dbsnp_vcf_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('allele_counts.tsv')
        ),
        kwargs={
            'chromosomes': 'autosomes'
        }
    )

    workflow.commandline(
        name='build_normal_wig',
        args=(
            'readCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(normal_bam_file),
            '>',
            mgd.TempOutputFile('normal.wig')
        )
    )

    workflow.commandline(
        name='build_tumour_wig',
        args=(
            'readCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(tumour_bam_file),
            '>',
            mgd.TempOutputFile('tumour.wig')
        )
    )

    workflow.commandline(
        name='build_gc_wig',
        args=(
            'gcCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(ref_genome_fasta_file),
            '>',
            mgd.TempOutputFile('gc.wig')
        )
    )

    workflow.commandline(
        name='build_mappability_wig',
        args=(
            'mapCounter',
            '-c', ','.join(chromosomes),
            mgd.InputFile(mappability_file),
            '>',
            mgd.TempOutputFile('mappability.wig')
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
            mgd.TempOutputFile('coverage.wig')
        ),
        kwargs={
            'target_file': exome_bed_file
        }
    )

    workflow.transform(
        name='run_titan',
        axes=('param_idx',),
        ctx={'mem': 8, 'mem_retry_increment': 4, 'num_retry': 3, 'threads': threads},
        func=tasks.run_titan,
        args=(
            mgd.TempInputFile('coverage.wig'),
            mgd.TempInputFile('allele_counts.tsv'),
            mgd.TempInputObj('init_params', 'param_idx'),
            mgd.TempOutputFile('run.tar.gz', 'param_idx'),
            mgd.TempSpace('titan_tmp', 'param_idx')
        ),
        kwargs={
            'is_exome': (exome_bed_file is not None),
            'sample': sample,
            'threads': threads
        }
    )

    workflow.transform(
        name='build_run_stats_file',
        func=tasks.build_run_stats_file,
        args=(
            mgd.TempInputFile('run.tar.gz', 'param_idx'),
            mgd.TempInputObj('init_params', 'param_idx'),
            mgd.TempOutputFile('stats.tsv')
        )
    )

    workflow.transform(
        name='build_output',
        func=tasks.build_final_results_file,
        args=(
            mgd.TempInputFile('coverage.wig'),
            mgd.TempInputFile('allele_counts.tsv'),
            mgd.TempInputFile('run.tar.gz', 'param_idx'),
            mgd.TempInputFile('stats.tsv'),
            mgd.OutputFile(out_file),
            mgd.TempSpace('build_results')
        )
    )

    return workflow


def create_allele_counts_workflow(
        normal_bam_file,
        tumour_bam_file,
        dbsnp_vcf_file,
        ref_genome_fasta_file,
        allele_counts_file,
        chromosomes='autosomes'):

    chromosomes = soil.utils.genome.load_bam_chromosome_lengths(normal_bam_file, chromosomes)

    sandbox = soil.utils.workflow.get_sandbox(['snpsift'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.subworkflow(
        name='call_snps',
        func=soil.wrappers.platypus.workflows.create_single_sample_workflow,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('normal.vcf.gz'),
        ),
        kwargs={
            'chromosomes': chromosomes,
            'split_size': int(1e7)
        }
    )

    workflow.commandline(
        name='annotate_dbsnp_status',
        ctx={'mem': 6, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'annotate',
            mgd.InputFile(dbsnp_vcf_file),
            mgd.TempInputFile('normal.vcf.gz'),
            '>',
            mgd.TempOutputFile('normal.dbsnp.vcf')
        )
    )

    workflow.commandline(
        name='annotate_variant_type',
        ctx={'mem': 6, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'varType',
            mgd.TempInputFile('normal.dbsnp.vcf'),
            '>',
            mgd.TempOutputFile('normal.dbsnp.vartype.vcf')
        )
    )

    workflow.commandline(
        name='filter_het_snps',
        ctx={'mem': 6, 'mem_retry_increment': 4, 'num_retry': 3},
        args=(
            'SnpSift',
            'filter',
            "isHet(GEN[0]) & ((exists ID) & ( ID =~ 'rs' )) & (exists SNP)",
            mgd.TempInputFile('normal.dbsnp.vartype.vcf'),
            '>',
            mgd.TempOutputFile('het.snps.vcf')
        )
    )

    workflow.transform(
        name='split_vcf',
        ctx={'mem': 6, 'mem_retry_increment': 4, 'num_retry': 3},
        func=tasks.split_vcf,
        args=(
            mgd.TempInputFile('het.snps.vcf'),
            mgd.TempOutputFile('split.vcf', 'split'),
            mgd.TempSpace('split_tmp')
        ),
        kwargs={
            'split_size': int(1e4)
        }
    )

    workflow.transform(
        name='get_allele_counts',
        axes=('split',),
        func=tasks.get_snv_allele_counts_for_vcf_targets,
        args=(
            mgd.InputFile(tumour_bam_file),
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('split.tsv', 'split')
        )
    )

    workflow.transform(
        name='merge_counts',
        func=tasks.merge_counts,
        args=(
            mgd.TempInputFile('split.tsv', 'split'),
            mgd.OutputFile(allele_counts_file)
        )
    )

    return workflow
