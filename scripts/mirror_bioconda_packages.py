from soil.utils.conda import packages

import subprocess

for package in packages.values():
    if package.channel != 'bioconda':
        continue

    pkg_spec = '{channel}/{name}/{version}'.format(**package.__dict__)

    cmd = [
        'anaconda',
        'copy',
        '--to-owner', 'soil',
        pkg_spec
    ]

    subprocess.call(cmd)
