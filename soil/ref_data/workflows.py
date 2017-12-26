import click
import os
import pypeliner
import pypeliner.managed as mgd
import yaml

import soil.ref_data.paths
import soil.ref_data.haplotype.workflows
import soil.utils.workflow
import soil.wrappers.bwa.tasks
import soil.wrappers.kallisto.tasks
import soil.wrappers.samtools.tasks
import soil.wrappers.star.tasks

import tasks


def create_ref_data_workflow(config, out_dir, cosmic=False, local_download=False, threads=1):
    """ Download and index reference data.
    """
    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='download',
        func=crete_download_ref_data_workflow,
        args=(
            config,
            out_dir
        ),
        kwargs={
            'cosmic': cosmic,
            'local_download': local_download
        }
    )

    workflow.subworkflow(
        name='index',
        func=create_index_ref_data_workflow,
        args=(
            config,
            out_dir
        ),
        kwargs={
            'cosmic': cosmic,
            'threads': threads
        }
    )

    return workflow


def crete_download_ref_data_workflow(config, out_dir, cosmic=False, local_download=False):
    """ Download reference files.

    This workflow mainly retrieves files from the internet. There are some light to moderately heavy computational tasks
    as well.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ref_data_paths = soil.ref_data.paths.SoilRefDataPaths(out_dir)

    with open(ref_data_paths.config_file, 'w') as fh:
        yaml.dump(config, fh)

    if cosmic:
        cosmic_user = click.prompt('Please enter COSMIC user ID')

        cosmic_password = click.prompt('Please enter COSMIC password', hide_input=True)

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    for key in config:
        if key.endswith('url') or key.endswith('urls'):
            workflow.setobj(obj=mgd.TempOutputObj(key), value=config[key])

    workflow.setobj(mgd.TempOutputObj('snpeff_url'), value=config['snpeff']['url'])

    workflow.subworkflow(
        name='download_ref_gene_annotations',
        func=_create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_gene_annotations_gtf_urls'),
            mgd.OutputFile(ref_data_paths.gene_annotations_gtf_file)
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.subworkflow(
        name='download_ref_genome',
        func=_create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_genome_fasta_urls'),
            mgd.TempOutputFile('raw_ref.fasta')
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.transform(
        name='lexsort_ref_genome',
        func=tasks.lex_sort_fasta,
        args=(
            mgd.TempInputFile('raw_ref.fasta'),
            mgd.OutputFile(ref_data_paths.genome_fasta_file)
        )
    )

    workflow.subworkflow(
        name='download_ref_proteome',
        func=_create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_proteome_fasta_urls'),
            mgd.TempOutputFile('raw_ref_prot.fasta')
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.transform(
        name='filter_bad_proteins',
        func=tasks.filter_bad_proiteins,
        args=(
            mgd.TempInputFile('raw_ref_prot.fasta'),
            mgd.OutputFile(ref_data_paths.proteome_fasta_file)
        )
    )

    workflow.subworkflow(
        name='download_ref_transcriptome',
        func=_create_download_decompress_concat_workflow,
        args=(
            mgd.TempInputObj('ref_transcriptome_fasta_urls'),
            mgd.OutputFile(ref_data_paths.transcriptome_fasta_file)
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.transform(
        name='download_dbsnp',
        ctx={'local': local_download},
        func=tasks.download,
        args=(
            mgd.TempInputObj('dbsnp_url'),
            mgd.OutputFile(ref_data_paths.dbsnp_vcf_file)
        )
    )

    if cosmic:
        workflow.subworkflow(
            name='download_cosmic',
            func=_create_download_cosmic_workflow,
            args=(
                config['cosmic']['ref_genome_version'],
                mgd.OutputFile(ref_data_paths.cosmic_vcf_file),
                cosmic_user,
                cosmic_password
            ),
            kwargs={
                'local_download': local_download
            }
        )

    workflow.transform(
        name='download_snpeff_db',
        ctx={'local': local_download},
        func=tasks.download,
        args=(
            mgd.TempInputObj('snpeff_url'),
            mgd.TempOutputFile('snpeff.zip')
        )
    )

    workflow.transform(
        name='unzip_snpeff',
        func=tasks.unzip_file,
        args=(
            mgd.TempInputFile('snpeff.zip'),
            mgd.OutputFile(os.path.join(os.path.dirname(ref_data_paths.snpeff_data_dir), 'done.txt')),
            mgd.TempSpace('snpeff_tmp')
        )
    )

    workflow.transform(
        name='download_genetic_map',
        ctx={'local': local_download},
        func=tasks.download,
        args=(
            mgd.TempInputObj('genetic_map_txt_url'),
            mgd.OutputFile(ref_data_paths.genetic_map_file)
        )
    )

    workflow.subworkflow(
        name='ref_haplotype_panel',
        func=soil.ref_data.haplotype.workflows.create_eagle_ref_data_workflow,
        args=(
            mgd.TempInputObj('haplotype_vcf_template_url'),
            mgd.OutputFile(ref_data_paths.haplotypes_bcf)
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.transform(
        name='download_iedb_mhc_one',
        ctx={'local': local_download},
        func=tasks.download,
        args=(
            mgd.TempInputObj('iedb_mhc_one_url'),
            mgd.TempOutputFile('mhc1.tar.gz')
        )
    )

    workflow.transform(
        name='extract_iedb_mhc_one',
        func=tasks.extract_tar_file,
        args=(
            mgd.TempInputFile('mhc1.tar.gz'),
            mgd.OutputFile(os.path.join(ref_data_paths.iedb_mhc_one_dir, 'extract.done'))
        )
    )

    workflow.transform(
        name='config_iedb_mhc_one',
        func=tasks.configure_iedb_module,
        args=(
            mgd.InputFile(os.path.join(ref_data_paths.iedb_mhc_one_dir, 'extract.done')),
            mgd.OutputFile(os.path.join(ref_data_paths.iedb_mhc_one_dir, 'configure.done'))
        )
    )

    workflow.transform(
        name='download_vep_cache',
        ctx={'local': local_download},
        func=tasks.download,
        args=(
            mgd.TempInputObj('vep_cache_url'),
            mgd.TempOutputFile('vep.tar.gz')
        )
    )

    workflow.transform(
        name='extract_vep_cache',
        func=tasks.extract_tar_file,
        args=(
            mgd.TempInputFile('vep.tar.gz'),
            mgd.OutputFile(os.path.join(ref_data_paths.vep_cache_dir, 'homo_sapiens', 'extract.done'))
        )
    )

    return workflow


def create_index_ref_data_workflow(out_dir, cosmic=False, threads=1):
    """ Create index files for references.

    This workflow is extremely compute and memory heavy. It should be run on a cluster with large memory nodes
    available.
    """
    ref_data_paths = soil.ref_data.paths.SoilRefDataPaths(out_dir)

    sandbox = soil.utils.workflow.get_sandbox(['bwa', 'bcftools', 'kallisto', 'samtools', 'star'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.commandline(
        name='link_bwa_ref',
        args=(
            'ln',
            mgd.InputFile(ref_data_paths.genome_fasta_file),
            mgd.OutputFile(ref_data_paths.bwa_genome_fasta_file)
        )
    )

    workflow.transform(
        name='bwa_index_ref_genome',
        ctx={'mem': 8, 'mem_retry_increment': 8, 'num_retry': 3},
        func=soil.wrappers.bwa.tasks.index,
        args=(
            mgd.InputFile(ref_data_paths.bwa_genome_fasta_file),
            mgd.OutputFile(ref_data_paths.bwa_genome_fasta_file + '.bwa_index.done')
        )
    )

    workflow.subworkflow(
        name='build_bwa_mappability_file',
        func=tasks.mappability_wrapper,
        args=(
            mgd.InputFile(ref_data_paths.bwa_genome_fasta_file + '.bwa_index.done'),
            mgd.OutputFile(ref_data_paths.genome_bwa_mappability_wig_file)
        ),
        kwargs={
            'k': 100,
            'max_map_qual': 60,
            'threads': threads
        }

    )

    workflow.commandline(
        name='link_star_ref',
        args=(
            'ln',
            mgd.InputFile(ref_data_paths.genome_fasta_file),
            mgd.OutputFile(ref_data_paths.star_genome_fasta_file)
        )
    )

    workflow.transform(
        name='star_index_ref_genome',
        ctx={'mem': 32, 'mem_retry_increment': 16, 'num_retry': 3, 'threads': threads},
        func=soil.wrappers.star.tasks.index,
        args=(
            mgd.InputFile(ref_data_paths.star_genome_fasta_file),
            mgd.InputFile(ref_data_paths.gene_annotations_gtf_file),
            mgd.OutputFile(ref_data_paths.star_genome_fasta_file + '.star_index.done')
        ),
        kwargs={
            'threads': threads
        }
    )

    workflow.transform(
        name='samtools_index_ref_genome',
        func=soil.wrappers.samtools.tasks.index_fasta,
        args=(
            mgd.InputFile(ref_data_paths.genome_fasta_file),
            mgd.OutputFile(ref_data_paths.genome_fasta_file + '.fai')
        )
    )

    workflow.transform(
        name='kallisto_index',
        ctx={'mem': 4, 'mem_retry_increment': 4, 'num_retry': 3},
        func=soil.wrappers.kallisto.tasks.build_index,
        args=(
            mgd.InputFile(ref_data_paths.transcriptome_fasta_file),
            mgd.OutputFile(ref_data_paths.kallisto_index_file)
        ),
        kwargs={
            'kmer_length': 31
        }
    )

    if cosmic:
        workflow.transform(
            name='index_cosmic',
            func=soil.wrappers.samtools.tasks.index_vcf,
            args=(
                mgd.InputFile(ref_data_paths.cosmic_vcf_file),
                mgd.OutputFile(ref_data_paths.cosmic_vcf_file + '.tbi')
            )
        )

    workflow.transform(
        name='index_dbsnp',
        func=soil.wrappers.samtools.tasks.index_vcf,
        args=(
            mgd.InputFile(ref_data_paths.dbsnp_vcf_file),
            mgd.OutputFile(ref_data_paths.dbsnp_vcf_file + '.tbi')
        )
    )

    return workflow


#=========================================================================
# Helper workflows
#=========================================================================


def _create_download_decompress_workflow(url, local_path, local_download=False):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(mgd.TempOutputObj('url'), value=url)

    workflow.transform(
        name='download',
        ctx={'local': local_download},
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


def _create_download_decompress_concat_workflow(urls, out_file, local_download=False):
    workflow = pypeliner.workflow.Workflow()

    local_files = []

    for i, url in enumerate(urls):
        local_files.append(mgd.TempFile('file_{}'.format(i)))

        workflow.setobj(mgd.TempOutputObj('url_{}'.format(i)), value=url)

        workflow.subworkflow(
            name='download_file_{}'.format(i),
            func=_create_download_decompress_workflow,
            args=(
                mgd.TempInputObj('url_{}'.format(i)),
                local_files[i].as_output(),
            ),
            kwargs={
                'local_download': local_download
            }
        )

    concat_args = ['cat', ] + [x.as_input() for x in local_files] + ['>', mgd.OutputFile(out_file)]

    workflow.commandline(
        name='concat',
        args=concat_args
    )

    return workflow


def _create_download_cosmic_workflow(
        ref_data_version,
        out_file,
        user,
        password,
        host='sftp-cancer.sanger.ac.uk',
        local_download=False):

    host_base_path = '/files/{}/cosmic/v83/VCF'.format(ref_data_version.lower())

    coding_host_path = '/'.join([host_base_path, 'CosmicCodingMuts.vcf.gz'])

    non_coding_host_path = '/'.join([host_base_path, 'CosmicNonCodingVariants.vcf.gz'])

    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(obj=mgd.TempOutputObj('coding_host_path'), value=coding_host_path)

    workflow.setobj(obj=mgd.TempOutputObj('non_coding_host_path'), value=non_coding_host_path)

    workflow.subworkflow(
        name='download_coding',
        func=_create_download_cosmic_file_subworkflow,
        args=(
            host,
            mgd.TempInputObj('coding_host_path'),
            user,
            password,
            mgd.TempOutputFile('coding.vcf.gz'),
        ),
        kwargs={
            'local_download': local_download
        }
    )

    workflow.subworkflow(
        name='download_non_coding',
        func=_create_download_cosmic_file_subworkflow,
        args=(
            host,
            mgd.TempInputObj('non_coding_host_path'),
            user,
            password,
            mgd.TempOutputFile('non_coding.vcf.gz'),
        ),
        kwargs={
            'local_download': local_download
        }
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


def _create_download_cosmic_file_subworkflow(host, host_path, user, password, out_file, local_download=False):
    sandbox = soil.utils.workflow.get_sandbox(['bcftools', 'samtools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.transform(
        name='download',
        ctx={'local': local_download},
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
        name='decompress',
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
