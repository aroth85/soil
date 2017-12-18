import soil.utils.package_data


def load_config_file(ref_genome_version):
    return soil.utils.package_data.load_data_file('ref_data/configs/{}.yaml'.format(ref_genome_version))
