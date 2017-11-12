import os
import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.star.workflows
import soil.wrappers.stringtie.tasks
import soil.wrappers.transdecoder.tasks


def create_assembly_from_fastq_workflow(
        fastq_file_1,
        fastq_file_2,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_bam_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=1):

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='align',
        func=soil.wrappers.star.workflows.create_align_workflow,
        args=(
            mgd.InputFile(fastq_file_1),
            mgd.InputFile(fastq_file_2),
            os.path.dirname(ref_genome_fasta_file),
            mgd.OutputFile(out_bam_file)
        ),
        kwargs={
            'add_xs_tag': True,
            'align_threads': threads,
            'sort_threads': threads
        }
    )

    workflow.subworkflow(
        name='assemble_transcriptome',
        func=create_assembly_from_bam_workflow,
        args=(
            mgd.InputFile(out_bam_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.InputFile(ref_gtf_file),
            mgd.OutputFile(out_cds_fasta_file),
            mgd.OutputFile(out_protein_fasta_file),
            mgd.OutputFile(out_gff_file),
            mgd.OutputFile(out_gtf_file),
        ),
        kwargs={
            'threads': threads,
        }
    )

    return workflow


def create_assembly_from_bam_workflow(
        bam_file,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['stringtie', 'transdecoder'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='run_stringtie',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        func=soil.wrappers.stringtie.tasks.assemble,
        args=(
            mgd.InputFile(bam_file),
            mgd.InputFile(ref_gtf_file),
            mgd.OutputFile(out_gtf_file),
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.transform(
        name='run_transdecoder',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=soil.wrappers.transdecoder.tasks.gtf_to_fasta,
        args=(
            mgd.InputFile(out_gtf_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(out_cds_fasta_file),
            mgd.OutputFile(out_protein_fasta_file),
            mgd.OutputFile(out_gff_file),
            mgd.TempSpace('transdecoder_tmp'),
        )
    )

    return workflow
