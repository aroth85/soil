import pkg_resources


def load_data_file(path):
    return pkg_resources.resource_filename(
        'soil',
        path
    )
