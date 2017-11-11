from setuptools import find_packages, setup

description = '''Wrapper library for bioinformatic tools.'''

setup(
    name='soil',
    version='0.1.0',
    description=description,
    author='Andrew Roth',
    author_email='andrewjlroth@gmail.com',
    url='https://github.com/aroth85/soil',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'soil-ref = soil.cli:ref',
            'soil-run = soil.cli:run',
        ]
    }
)
