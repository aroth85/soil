import os


def find(name, path):
    """ Find a file in a given path. Useful for locating executables in the conda environment.

    :param name: Name of file to find.
    :param path: Path to search for file in.
    :returns: Absolute path to the file.
    """
    for root, _, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
