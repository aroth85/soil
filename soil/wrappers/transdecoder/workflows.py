import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow

import tasks


def create_transdecoder_workflow(
        in_gtf_file,
        ref_gtf_file,
        ref_genome_fasta_file,
        out_alignment_gff_file,
        out_cdna_fasta_file,
        out_cds_fasta_file,
        out_protein_fasta_file):

    sandbox = soil.utils.workflow.get_sandbox(['transdecoder', ])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='convert_gtf_to_cdna_fasta',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.convert_gtf_to_cdna_fasta,
        args=(
            mgd.InputFile(in_gtf_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(out_cdna_fasta_file),
        )
    )

    workflow.transform(
        name='convert_gtf_to_gff',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.convert_gtf_to_gff_file,
        args=(
            mgd.InputFile(in_gtf_file),
            mgd.TempOutputFile('ref.gff'),
        )
    )

    workflow.transform(
        name='run_transdecoder',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.run_transdecoder,
        args=(
            mgd.InputFile(out_cdna_fasta_file),
            mgd.OutputFile(out_cds_fasta_file),
            mgd.OutputFile(out_protein_fasta_file),
            mgd.TempOutputFile('transdecoder.gff'),
            mgd.TempSpace('transdecoder_tmp'),
        ),

    )

    workflow.transform(
        name='buil_alignment_gff',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=tasks.build_alignment_gff,
        args=(
            mgd.InputFile(out_cdna_fasta_file),
            mgd.TempInputFile('transdecoder.gff'),
            mgd.TempInputFile('ref.gff'),
            mgd.OutputFile(out_alignment_gff_file),
        )
    )

    return workflow
