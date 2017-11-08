import click

import soil.wrappers.bwa.cli
import soil.wrappers.mixcr.cli
import soil.wrappers.platypus.cli
import soil.wrappers.strelka.cli
import soil.wrappers.varscan.cli


@click.group()
def run():
    pass

run.add_command(soil.wrappers.bwa.cli.bwa)
run.add_command(soil.wrappers.mixcr.cli.mixcr)
run.add_command(soil.wrappers.platypus.cli.platypus)
run.add_command(soil.wrappers.strelka.cli.strelka)
run.add_command(soil.wrappers.varscan.cli.varscan)


@click.group()
def pipeline():
    pass
