import click

import soil.pipelines.dna_db.cli
import soil.pipelines.rna_assembly.cli

import soil.ref_data.cli

import soil.wrappers.bwa.cli
import soil.wrappers.eagle.cli
import soil.wrappers.mixcr.cli
import soil.wrappers.msgf_plus.cli
import soil.wrappers.optitype.cli
import soil.wrappers.platypus.cli
import soil.wrappers.star.cli
import soil.wrappers.strelka.cli
import soil.wrappers.topiary.cli
import soil.wrappers.titan.cli
import soil.wrappers.transdecoder.cli


@click.group()
def pipeline():
    pass


pipeline.add_command(soil.pipelines.rna_assembly.cli.rna_assembly)
pipeline.add_command(soil.pipelines.dna_db.cli.dna_db)


@click.group()
def ref():
    """ Tools for handling reference data files used by soil.
    """
    pass


ref.add_command(soil.ref_data.cli.download)
ref.add_command(soil.ref_data.cli.index)
ref.add_command(soil.ref_data.cli.mappability)
ref.add_command(soil.ref_data.cli.show_config)


@click.group()
def run():
    pass


run.add_command(soil.wrappers.bwa.cli.bwa)
run.add_command(soil.wrappers.eagle.cli.eagle)
run.add_command(soil.wrappers.mixcr.cli.mixcr)
run.add_command(soil.wrappers.msgf_plus.cli.msgf_plus)
run.add_command(soil.wrappers.optitype.cli.optitype)
run.add_command(soil.wrappers.platypus.cli.platypus)
run.add_command(soil.wrappers.star.cli.star)
run.add_command(soil.wrappers.strelka.cli.strelka)
run.add_command(soil.wrappers.topiary.cli.topiary)
run.add_command(soil.wrappers.titan.cli.titan)
run.add_command(soil.wrappers.transdecoder.cli.transdecoder)
