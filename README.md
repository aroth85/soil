# Southampton Oxford Immunomics Library (SOIL)

A library to support colloborative immunomics research in Southampton and Oxford.
SOIL provides collection of wrappers and workflows to support analysis across multiple types of omics data.

# Installation

The SOIL package is currently under development so the installation process is currently a bit convoluted.
In future this will be simplified using conda packaging.

You will need to have `conda` installed. 
If you don not the best way to get conda is to use the [Miniconda](https://conda.io/miniconda.html) installer.

Once you have `conda` installed it is good practice to create a `conda` environment for SOIL with Python 2.

```
conda create -n soil python=2.7
source activate soil
```

Next you will need to get [SOIL](https://github.com/aroth85/pypeliner) and an experimental fork of [Pypeliner](https://github.com/aroth85/pypeliner) using git.
The following assumes you are in a directory where you want to put the code.
First install Pypeliner.

```
git clone https://github.com/aroth85/pypeliner.git
cd pypeliner
python setup.py develop
conda install --yes --file requirements.txt
```

Next install SOIL, again assuming you are in the directory where the code resides i.e. `cd ..` from the previous steps.

```
git clone git@github.com:aroth85/soil.git
cd soil
python setup.py develop
conda install --yes --file requirements.txt -c bioconda
pip install pyopenms
```

If all worked then you should have `soil-run` on the path and it won't raise errors.

> Note: If you want to run on a cluster install the `drmaa` package using conda.

# Using SOIL wrappers

The simples way to use soil is through the `soil-run` command which provides wrappers for several useful programs.
Wrappers do several things:

1. Perform the required steps to get a useful output from the analysis. 
For example the `soil-run bwa mem` wrapper runs `bwa`, sorts the output and marks duplicates.

2. Install the required software on the fly.
SOIL uses the concept of sandboxes for tasks or workflows to ensure the correct software is available.
The software is installed in a `conda` environment on the fly.

3. Supports execution on a local machine or a cluster.


Currentyl there is support for
- bwa
- mixcr
- msgf_plus
- optitype
- platypus
- star
- strelka
- transdecoder

To see the help for a wrapper use the `--help` flag.
For example

```
soil-run optitype --help
```

would show the following

```
Usage: soil-run optitype [OPTIONS]

Options:
  -b, --bam-file PATH          Path to BAM file to analyze. Should be a normal (non-malignant) sample if available.
                               [required]
  -o, --out-file PATH          Path where output with HLA information will be written.  [required]
  --rna                        Set this flag if the data is from RNA (WTS). Otherwise assumed to be DNA (WGS/WES).
  -th, --threads INTEGER       Number of threads used in parallel steps of workflow.
  -wd, --working-dir PATH      Working directory for runner. Analysis will fail if it exists unless --resume flag is
                               used.Will be deleted when analysis is finished.  [required]
  -mj, --max-jobs INTEGER      Maximum number of jobs to run.
  -ns, --native-spec TEXT      String specifying cluster submission parameters. Special values are {mem} for memory
                               requests and {threads} for thread requests. This can be set globally by through the
                               SOIL_NATIVE_SPEC environment variable See online help for examples.
  -sb, --submit [drmaa|local]  Job submission strategy. Use local to run on host machine or drmaa to submit to a
                               cluster.
  --resume                     Set this flag if an analysis was interrupted and you would like to resume. Only has an
                               effect if the working directory exists.
  --help                       Show this message and exit.
```

Some programs may have subcommands, for example `soil-run bwa mem`.

## Common options

All wrappers will have the following options

```
  -wd, --working-dir PATH      Working directory for runner. Analysis will fail if it exists unless --resume flag is
                               used.Will be deleted when analysis is finished.  [required]
  -mj, --max-jobs INTEGER      Maximum number of jobs to run.
  -ns, --native-spec TEXT      String specifying cluster submission parameters. Special values are {mem} for memory
                               requests and {threads} for thread requests. This can be set globally by through the
                               SOIL_NATIVE_SPEC environment variable See online help for examples.
  -sb, --submit [drmaa|local]  Job submission strategy. Use local to run on host machine or drmaa to submit to a
                               cluster.
  --resume                     Set this flag if an analysis was interrupted and you would like to resume. Only has an
                               effect if the working directory exists.
  --help                       Show this message and exit.
```

### Working directory -wd

Under the hood SOIL uses Pypeliner to run workflows. 
These workflows require a working directory to write temporary files and logs.
This flag controls where that directory is stored.
Some important points

1. The directory should be on a disk with enough space for the analysis.
2. The directory should be accessible to compute nodes on the cluster.
3. The directory is automatically cleaned up after analysis finishes successfully.
4. If the analysis fails, for example because the cluster goes down the `--resume` flag can be used to restart from where things left off.

## Max jobs -mj

SOIL will try to parallelize jobs where possible.
This flag controls the maximum number of concurrent jobs which will be run.
When running locally this should be less equal the number of cores available.
For clusters it is less important as the queue system will control things if to many jobs get launched.

> Note: If you are running locally and using a wrapper with a thread command then you may overload the system.
> SOIL is not aware of how many threads a job is using when it computes max jobs, so you can request more threads than jobs.
> This is not an issue on clusters as the submit system takes care of the issues.

##  Native spec -ns

This flag is only relevant if the wrapper is being run on a cluster.
This is equivalent to the arguments that would be put in a qsub script.
For example

```
-ns "-cwd -V -pe smp {threads} -q byslot.q -l mem_free={mem}G,h_vmem={mem}G"
```

requests the analysis be launched in the current working directory (`-cwd`) inheriting the current environment (`-V`) etc.
Of note are the special strings `{threads}` and `{mem}` which control the number of threads requested and memory.
These need to be setup approriately by the user for the cluster.

Note that it is possible to set the environment variable `SOIL_NATIVE_SPEC`, in which case the `-ns` flag is not needed.

## Submit -sb

SOIL can run either on a local machine or on a cluster (currently only Grid Engine is tested).
The `-sb` flag chooses which, with `local` running locally and `drmaa` on a cluster via the DRMAA API.

> Note: If you want to run on a cluster install the `drmaa` package using conda.