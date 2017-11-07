"""
Created on 7 Nov 2017

@author: Andrew Roth
"""
import click

import soil.wrappers.strelka.cli
import soil.wrappers.varscan.cli


@click.group()
def run():
    pass

run.add_command(soil.wrappers.strelka.cli.strelka)
run.add_command(soil.wrappers.varscan.cli.varscan)


@click.group()
def pipeline():
    pass
