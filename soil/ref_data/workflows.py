import click
import os
import pypeliner
import pypeliner.managed as mgd
import yaml

import soil.utils.workflow
import soil.wrappers.bwa.tasks
import soil.wrappers.samtools.tasks

import tasks


def create_ref_data_workflow(ref_genome_version, out_dir, cosmic=False):
    if ref_genome_version == 'GRCh37':
        import soil.ref_data.configs.GRCh37

        config = soil.ref_data.configs.GRCh37.get_config()

    else:
        raise Exception('{} is not a supported reference'.format(ref_genome_version))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, 'config.yaml'), 'w') as fh:
        yaml.dump(config, fh)

    if cosmic:
        cosmic_user = click.prompt('Please enter COSMIC user ID')

        cosmic_password = click.prompt('Please enter COSMIC password', hide_input=True)

        cosmic_file = os.path.join(out_dir, 'cosmic.vcf.gz')

    dbsnp_file = os.path.join(out_dir, 'dbsnp.vcf.gz')

    ref_gene_annotations_gtf_file = os.path.join(out_dir, 'ref_gene_annotations.gtf')

    ref_genome_fasta_file = os.path.join(out_dir, 'ref_genome.fasta')

    ref_proteome_fasta_file = os.path.join(out_dir, 'ref_proteome.fasta')

    sandbox = soil.utils.workflow.get_sandbox(['bwa', 'bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.subworkflow(
        name='download_ref_gene_annotations',
        func=create_download_decompress_workflow,
        args=(
            config['ref_gene_annotations_gtf_url'],
            mgd.OutputFile(ref_gene_annotations_gtf_file)
        )
    )

    workflow.subworkflow(
        name='download_ref_genome',
        func=create_download_decompress_workflow,
        args=(
            config['ref_genome_fasta_url'],
            mgd.TempOutputFile('raw_ref.fasta'),
        )
    )

    workflow.transform(
        name='lexsort_ref_genome',
        func=tasks.lex_sort_fasta,
        args=(
            mgd.TempInputFile('raw_ref.fasta'),
            mgd.OutputFile(ref_genome_fasta_file),
        )
    )
    workflow.transform(
        name='bwa_index_ref_genome',
        func=soil.wrappers.bwa.tasks.index,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(ref_genome_fasta_file + '.bwa_index.done'),
        )
    )

    workflow.transform(
        name='samtools_index_ref_genome',
        func=soil.wrappers.samtools.tasks.index_fasta,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.OutputFile(ref_genome_fasta_file + '.fai'),
        )
    )

    workflow.subworkflow(
        name='download_ref_proteome',
        func=create_download_decompress_workflow,
        args=(
            config['ref_proteome_fasta_url'],
            mgd.OutputFile(ref_proteome_fasta_file)
        )
    )

    workflow.subworkflow(
        name='download_dbsnp',
        func=create_download_workflow,
        args=(
            config['dbsnp_url'],
            mgd.OutputFile(dbsnp_file)
        )
    )

    workflow.transform(
        name='index_dbsnp',
        func=soil.wrappers.samtools.tasks.index_vcf,
        args=(
            mgd.InputFile(dbsnp_file),
            mgd.OutputFile(dbsnp_file + '.tbi'),
        )
    )

    if cosmic:
        workflow.subworkflow(
            name='download_cosmic',
            func=create_download_cosmic_workflow,
            args=(
                ref_genome_version,
                mgd.OutputFile(cosmic_file),
                cosmic_user,
                cosmic_password,
            )
        )

    return workflow


def create_download_workflow(url, local_path):
    """ Download a file from a url.

    Encapsulate the simple task in a worklfow to avoid pypeliner rerunning.

    :param url: URL of file to download.
    :param local_path: Path where downloaded file will be placed.
    """
    workflow = pypeliner.workflow.Workflow()

    if not os.path.exists(local_path):
        workflow.transform(
            name='download',
            func=tasks.download,
            args=(
                url,
                mgd.OutputFile(local_path),
            ),
        )

    return workflow


def create_download_decompress_workflow(url, local_path):
    workflow = pypeliner.workflow.Workflow()

    if not os.path.exists(local_path):
        workflow.transform(
            name='download',
            func=tasks.download,
            args=(
                url,
                mgd.TempOutputFile('download.gz'),
            ),
        )

        workflow.transform(
            name='decompress',
            func=tasks.decompress,
            args=(
                mgd.TempInputFile('download.gz'),
                local_path
            )
        )

    return workflow


def create_download_cosmic_workflow(ref_data_version, out_file, user, password, host='sftp-cancer.sanger.ac.uk'):
    host_base_path = '/files/{}/cosmic/v83/VCF'.format(ref_data_version.lower())

    if not os.path.exists(out_file):
        sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

        workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

        workflow.subworkflow(
            name='download_coding',
            func=create_download_cosmic_file_subworkflow,
            args=(
                host,
                '/'.join([host_base_path, 'CosmicCodingMuts.vcf.gz']),
                user,
                password,
                mgd.TempOutputFile('coding.vcf.gz'),
            )
        )

        workflow.subworkflow(
            name='download_non_coding',
            func=create_download_cosmic_file_subworkflow,
            args=(
                host,
                '/'.join([host_base_path, 'CosmicNonCodingVariants.vcf.gz']),
                user,
                password,
                mgd.TempOutputFile('non_coding.vcf.gz'),
            )
        )

        workflow.transform(
            name='merge_files',
            func=soil.wrappers.samtools.tasks.concatenate_vcf,
            args=(
                [mgd.TempInputFile('coding.vcf.gz'), mgd.TempInputFile('non_coding.vcf.gz')],
                mgd.OutputFile(out_file)
            ),
            kwargs={
                'allow_overlap': True,
                'index_file': mgd.OutputFile(out_file + '.tbi')
            }
        )

    else:
        workflow = pypeliner.workflow.Workflow()

    return workflow


def create_download_cosmic_file_subworkflow(host, host_path, user, password, out_file):
    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='download_non_coding',
        func=tasks.download_from_sftp,
        args=(
            host,
            host_path,
            mgd.TempOutputFile('file.vcf.gz'),
            user,
            password
        )
    )

    workflow.transform(
        name='decompress_coding',
        func=tasks.decompress,
        args=(
            mgd.TempInputFile('file.vcf.gz'),
            mgd.TempOutputFile('file.vcf')
        )
    )

    workflow.transform(
        name='bgzip',
        func=soil.wrappers.samtools.tasks.compress_vcf,
        args=(
            mgd.TempInputFile('file.vcf'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow
