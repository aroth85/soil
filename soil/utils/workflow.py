def flatten_input(files):
    """ Takes a collection and returns a list. Useful for writing functions that may be used as one-to-one or
    many-to-one actions by pypeliner.

    :param files: A collection, usually of file paths.
    :returns: A flattened list.

    >>> flatten_input({'a': '/foo/bar/a', 'b': '/foo/bar/b'})
    ['/foo/bar/a', '/foo/bar/b']

    >>> flatten_input('/foo/bar/a')
    [''/foo/bar/a',]

    """
    if type(files) == dict:
        parsed_files = [files[x] for x in sorted(files)]

    elif type(files) == str:
        parsed_files = [files, ]

    else:
        parsed_files = []

        for x in files:
            if type(x) == dict:
                parsed_files.extend([x[y] for y in sorted(x)])
            else:
                parsed_files.append(x)

    return parsed_files