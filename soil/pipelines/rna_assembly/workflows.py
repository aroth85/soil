import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.stringtie.tasks
import soil.wrappers.transdecoder.workflows


def create_assembly_workflow(
        bam_file,
        ref_genome_fasta_file,
        ref_gtf_file,
        out_cds_fasta_file,
        out_protein_fasta_file,
        out_gff_file,
        out_gtf_file,
        threads=1):

    sandbox = soil.utils.workflow.get_sandbox(['stringtie'])

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

    workflow.subworkflow(
        name='run_transdecoder',
        func=soil.wrappers.transdecoder.workflows.create_transdecoder_workflow,
        args=(
            mgd.InputFile(out_gtf_file),
            mgd.InputFile(ref_gtf_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(out_gff_file),
            mgd.TempOutputFile('cdna.fasta'),
            mgd.OutputFile(out_cds_fasta_file),
            mgd.OutputFile(out_protein_fasta_file),
        )
    )

    return workflow
