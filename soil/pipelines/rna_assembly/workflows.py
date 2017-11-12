import pypeliner
import pypeliner.managed as mgd

import soil.utils.workflow
import soil.wrappers.stringtie.tasks
import soil.wrappers.transdecoder.tasks


def create_assembly_workflow(bam_file, ref_genome_fasta_file, ref_gtf_file, out_prefix, threads=1):
    cds_fasta_file = out_prefix + '.cds.fasta'

    protein_fasta_file = out_prefix + '.prot.fasta'

    gff_file = out_prefix + '.gff3'

    stringtie_gtf_file = out_prefix + '.gtf'

    sandbox = soil.utils.workflow.get_sandbox(['stringtie', 'transdecoder'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='run_stringtie',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3, 'threads': threads},
        func=soil.wrappers.stringtie.tasks.assemble,
        args=(
            mgd.InputFile(bam_file),
            mgd.InputFile(ref_gtf_file),
            mgd.OutputFile(stringtie_gtf_file),
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
            mgd.InputFile(stringtie_gtf_file),
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(cds_fasta_file),
            mgd.OutputFile(protein_fasta_file),
            mgd.OutputFile(gff_file),
            mgd.TempSpace('transdecoder_tmp'),
        )
    )

    return workflow
