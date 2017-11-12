import click
import os
import pypeliner
import pypeliner.managed as mgd
import yaml

import soil.utils.workflow
import soil.wrappers.bwa.tasks
import soil.wrappers.samtools.tasks
import soil.wrappers.star.tasks

import tasks


def create_ref_data_workflow(ref_genome_version, out_dir, cosmic=False, threads=1):
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

    ref_genome_fasta_file = os.path.join(out_dir, 'ref_genome.fa')

    ref_proteome_fasta_file = os.path.join(out_dir, 'ref_proteome.fa')

    sandbox = soil.utils.workflow.get_sandbox(['bwa', 'bcftools', 'samtools', 'star'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    for key in config:
        if key.endswith('url') or key.endswith('urls'):
            workflow.setobj(obj=mgd.TempOutputObj(key), value=config[key])

    workflow.subworkflow(
        name='download_ref_gene_annotations',
        func=create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_gene_annotations_gtf_urls'),
            mgd.OutputFile(ref_gene_annotations_gtf_file)
        )
    )

    workflow.subworkflow(
        name='download_ref_genome',
        func=create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_genome_fasta_urls'),
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
        name='star_index_ref_genome',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': threads},
        func=soil.wrappers.star.tasks.index,
        args=(
            mgd.InputFile(ref_genome_fasta_file),
            mgd.InputFile(ref_gene_annotations_gtf_file),
            mgd.OutputFile(ref_genome_fasta_file + '.star_index.done'),
        ),
        kwargs={
            'threads': threads,
        }
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
        func=create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_proteome_fasta_urls'),
            mgd.TempOutputFile('raw_ref_prot.fasta')
        )
    )

    workflow.transform(
        name='filter_bad_proteins',
        func=tasks.filter_bad_proiteins,
        args=(
            mgd.TempInputFile('raw_ref_prot.fasta'),
            mgd.OutputFile(ref_proteome_fasta_file)
        )
    )

    workflow.subworkflow(
        name='download_dbsnp',
        func=create_download_workflow,
        args=(
            mgd.TempInputObj('dbsnp_url'),
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

    workflow.setobj(mgd.TempOutputObj('url'), value=url)

    workflow.transform(
        name='download',
        func=tasks.download,
        args=(
            mgd.TempInputObj('url'),
            mgd.OutputFile(local_path),
        ),
    )

    return workflow


def create_download_decompress_workflow(url, local_path):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(mgd.TempOutputObj('url'), value=url)

    workflow.transform(
        name='download',
        func=tasks.download,
        args=(
            mgd.TempInputObj('url'),
            mgd.TempOutputFile('download'),
        ),
    )

    workflow.transform(
        name='decompress',
        func=tasks.decompress,
        args=(
            mgd.TempInputFile('download'),
            mgd.OutputFile(local_path),
        )
    )

    return workflow


def create_download_decompress_concat_workflow(urls, out_file):
    workflow = pypeliner.workflow.Workflow()

    if os.path.exists(out_file):
        return workflow

    local_files = []

    for i, url in enumerate(urls):
        local_files.append(mgd.TempFile('file_{}'.format(i)))

        workflow.setobj(mgd.TempOutputObj('url_{}'.format(i)), value=url)

        workflow.subworkflow(
            name='download_file_{}'.format(i),
            func=create_download_decompress_workflow,
            args=(
                mgd.TempInputObj('url_{}'.format(i)),
                local_files[i].as_output(),
            )
        )

    concat_args = ['cat', ] + [x.as_input() for x in local_files] + ['>', mgd.OutputFile(out_file)]

    workflow.commandline(
        name='concat_fastas',
        args=concat_args
    )

    return workflow


def create_download_cosmic_workflow(ref_data_version, out_file, user, password, host='sftp-cancer.sanger.ac.uk'):
    host_base_path = '/files/{}/cosmic/v83/VCF'.format(ref_data_version.lower())

    coding_host_path = '/'.join([host_base_path, 'CosmicCodingMuts.vcf.gz'])

    non_coding_host_path = '/'.join([host_base_path, 'CosmicNonCodingVariants.vcf.gz'])

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(obj=mgd.TempOutputObj('coding_host_path'), value=coding_host_path)

    workflow.setobj(obj=mgd.TempOutputObj('non_coding_host_path'), value=non_coding_host_path)

    workflow.subworkflow(
        name='download_coding',
        func=create_download_cosmic_file_subworkflow,
        args=(
            host,
            mgd.TempInputObj('coding_host_path'),
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
            mgd.TempInputObj('non_coding_host_path'),
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


def create_ref_genome_workflow(urls, out_file):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(obj=mgd.TempOutputObj('urls'), value=urls)

    workflow.subworkflow(
        name='download_ref_fasta_files',
        func=create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('urls'),
            mgd.TempOutputFile('raw_ref.fasta')
        )
    )

    workflow.transform(
        name='lexsort_ref_genome',
        func=tasks.lex_sort_fasta,
        args=(
            mgd.TempInputFile('raw_ref.fasta'),
            mgd.OutputFile(out_file),
        )
    )

    return workflow
