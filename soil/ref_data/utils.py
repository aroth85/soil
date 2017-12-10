import pkg_resources


def load_config_file(ref_genome_version):
    return pkg_resources.resource_filename(
        'soil',
        'ref_data/configs/{}.yaml'.format(ref_genome_version)
    )
