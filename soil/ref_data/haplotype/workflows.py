import pandas as pd
import pypeliner
import pypeliner.managed as mgd

import soil.ref_data.tasks
import soil.utils.package_data
import soil.utils.workflow
import soil.wrappers.bcftools.tasks

import tasks


def create_eagle_ref_data_workflow(vcf_url_template, out_file):
    chrom_map_file = soil.utils.package_data.load_data_file('ref_data/data/GRCh37/chrom_map.tsv')

    chrom_map = pd.read_csv(chrom_map_file, sep='\t')

    chrom_map = chrom_map[chrom_map['ncbi'] != 'MT']

    chrom_map['url'] = chrom_map['ncbi'].apply(lambda x: vcf_url_template.format(chrom=x))

    vcf_urls = chrom_map['url'].to_dict()

    sandbox = soil.utils.workflow.get_sandbox(['bcftools'])

    workflow = pypeliner.workflow.Workflow(default_sandbox=sandbox)

    workflow.setobj(obj=mgd.TempOutputObj('vcf_url', 'chrom'), value=vcf_urls)

    workflow.transform(
        name='download_vcf_files',
        axes=('chrom',),
        func=soil.ref_data.tasks.download,
        args=(
            mgd.TempInputObj('vcf_url', 'chrom'),
            mgd.TempOutputFile('raw.vcf.gz', 'chrom')
        )
    )

    workflow.transform(
        name='write_chrom_map',
        func=tasks.write_chrom_map_file,
        args=(
            mgd.InputFile(chrom_map_file),
            mgd.TempOutputFile('chrom_map.tsv')
        )
    )

    workflow.transform(
        name='rename_chroms',
        axes=('chrom',),
        func=soil.wrappers.bcftools.tasks.rename_chroms,
        args=(
            mgd.TempInputFile('chrom_map.tsv'),
            mgd.TempInputFile('raw.vcf.gz', 'chrom'),
            mgd.TempOutputFile('renamed.bcf', 'chrom')
        )
    )

    workflow.transform(
        name='concat_vcfs',
        func=soil.wrappers.bcftools.tasks.concatenate_vcf,
        args=(
            mgd.TempInputFile('renamed.bcf', 'chrom'),
            mgd.OutputFile(out_file)
        )
    )

    workflow.commandline(
        name='index',
        args=(
            'bcftools',
            'index',
            mgd.InputFile(out_file),
            '-o', mgd.OutputFile(out_file + '.csi')
        )
    )

    return workflow
