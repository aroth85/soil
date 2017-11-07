"""
Created on 7 Nov 2017

@author: Andrew Roth
"""
import click

import soil.wrappers.strelka.cli


@click.group()
def run():
    pass

run.add_command(soil.wrappers.strelka.cli.strelka)


@click.group()
def pipeline():
    pass
