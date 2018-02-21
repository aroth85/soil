import pypeliner
import pypeliner.managed as mgd
import pypeliner.sandbox

import soil.utils.workflow
import soil.wrappers.bwa.workflows
import soil.wrappers.sambamba.tasks
import soil.wrappers.strelka.workflows

import tasks


def create_custom_dna_proteome_from_fastq_workflow(
        normal_fastq_file_1,
        normal_fastq_file_2,
        tumour_fastq_file_1,
        tumour_fastq_file_2,
        ref_genome_fasta_file,
        ref_proteome_fasta_file,
        normal_bam_file,
        tumour_bam_file,
        custom_proteome_file,
        strelka_file,
        genome_version='GRCh37',
        is_exome=False,
        pyensembl_cache_dir=None,
        threads=1):

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='align_normal',
        func=create_align_workflow,
        args=(
            mgd.InputFile(normal_fastq_file_1),
            mgd.InputFile(normal_fastq_file_2),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(normal_bam_file)
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.subworkflow(
        name='align_tumour',
        func=create_align_workflow,
        args=(
            mgd.InputFile(tumour_fastq_file_1),
            mgd.InputFile(tumour_fastq_file_2),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(tumour_bam_file)
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.subworkflow(
        name='create_db',
        func=create_custom_proteom_from_bam_workflow,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.InputFile(ref_proteome_fasta_file),
            mgd.OutputFile(custom_proteome_file),
            mgd.OutputFile(strelka_file)
        ),
        kwargs={
            'genome_version': genome_version,
            'is_exome': is_exome,
            'pyensembl_cache_dir': pyensembl_cache_dir
        }
    )

    return workflow


def create_custom_proteom_from_bam_workflow(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        ref_proteome_fasta_file,
        custom_proteome_file,
        strelka_file,
        genome_version='GRCh37',
        is_exome=False,
        pyensembl_cache_dir=None,):

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.subworkflow(
        name='call_somatic_ssms',
        func=soil.wrappers.strelka.workflows.create_somatic_workflow,
        args=(
            mgd.InputFile(normal_bam_file),
            mgd.InputFile(tumour_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(strelka_file)
        ),
        kwargs={
            'is_exome': is_exome
        }
    )

    workflow.subworkflow(
        name='create_db',
        func=create_db_workflow,
        args=(
            mgd.InputFile(strelka_file),
            mgd.InputFile(ref_proteome_fasta_file),
            mgd.OutputFile(custom_proteome_file),
        ),
        kwargs={
            'genome_version': genome_version,
            'pyensembl_cache_dir': pyensembl_cache_dir
        }
    )

    return workflow


def create_align_workflow(fastq_file_1, fastq_file_2, ref_genome_fasta_file, out_bam_file, threads=1):
    sandbox = soil.utils.workflow.get_sandbox(['sambamba', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.subworkflow(
        name='align',
        func=soil.wrappers.bwa.workflows.create_align_workflow,
        args=(
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.TempOutputFile('aligned.bam')
        ),
        kwargs={
            'align_threads': threads,
            'sort_threads': threads
        }
    )

    workflow.transform(
        name='mark_dups',
        func=soil.wrappers.sambamba.tasks.markdups,
        args=(
            mgd.TempInputFile('aligned.bam'),
            mgd.OutputFile(out_bam_file),
            mgd.TempSpace('mark_dups_tmp')
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.commandline(
        name='index',
        args=(
            'samtools', 'index',
            mgd.InputFile(out_bam_file),
            mgd.OutputFile(out_bam_file + '.bai'),
        )
    )

    return workflow


def create_db_workflow(
        in_file,
        ref_proteome_fasta_file,
        out_file,
        genome_version='GRCh37',
        pyensembl_cache_dir=None):

    sandbox = pypeliner.sandbox.CondaSandbox(pip_packages=['varcode'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='clean_ref_fasta',
        func=tasks.clean_ref_proteome_ids,
        args=(
            mgd.InputFile(ref_proteome_fasta_file),
            mgd.TempOutputFile('ref.fasta')
        )
    )

    workflow.transform(
        name='build_variant_table',
        func=tasks.build_variant_table,
        args=(
            mgd.InputFile(in_file),
            mgd.TempOutputFile('variant_table.tsv.gz')
        ),
        kwargs={
            'genome_version': genome_version,
            'pyensembl_cache_dir': pyensembl_cache_dir
        }
    )

    workflow.transform(
        name='build_variant_fasta',
        func=tasks.build_variant_fasta,
        args=(
            mgd.TempInputFile('variant_table.tsv.gz'),
            mgd.TempOutputFile('var.fasta')
        )
    )

    workflow.commandline(
        name='build_db',
        args=(
            'cat',
            mgd.TempInputFile('ref.fasta'),
            mgd.TempInputFile('var.fasta'),
            '>',
            mgd.OutputFile(out_file)
        )
    )

    return workflow
