# TARGETSEQ
A nanopore pipeline based on the artic-ncov2019 pipeline

## Requirements
This pipeline wrapper requires python v3.6 or higher (Pyhton 3.7 is preferred)

Several dependencies are noted in the requirements.yml file and are installed automatically when creating the conda environment.
A specific commit of Jvarkit was taken from (http://lindenb.github.io/jvarkit/SAM4WebLogo.html, https://github.com/lindenb/jvarkit.git) 
and permanently added to the repo

This pipeline also requires a "sample_names.csv" to indicate which barcodes correspond to which samples for the demultiplexing step.
A template of this file is included in the repo.

This pipeline assumes that the primer scheme bed files are of a specific format:
 tab seperated ".bed" file 
 
with the column headings: "genome", "start", "end", "Primer_ID" and "number"
start is the start position of the primer, relative to the reference genome, using a zero based index (pos 1 = 0)
Primer names need to include "LEFT" or "RIGHT", using "_" as a delimeter to refer to Fwd and Rev primers

If you need to use a primer scheme that is not included here, please create an issue on this github page and I will add it to the repo for you

This pipeline will call the artic-ncov2019 pipeline which needs to be cloned seperately with its own environment

## Step 1
Download and install the 64-bit Python 3.7 version of Miniconda/Anaconda

## Step 2 clone the pipeline repo
`git clone --recursive https://github.com/PhilliVanilli/target.git`

 change into the repo directory
 
 `cd target`

## Step 3 create the conda environment
`conda env create -f requirements_v5.yml`

if pip ssl issues, create a HOME/.pip/pip.conf file with following text

[global]
trusted-host = pypi.python.org
               pypi.org
               files.pythonhosted.org

# Step 3.1 activate the conda env
activate the environment

`conda activate target`

## Step 3.2
install sam4weblogo from jvarkit folder, this tool converts a bam/sam file to a multi sequence alignment

` cd jvarkit`
`./gradlew sam4weblogo`

go back to base environment and root folder
`cd ..`
`cd ..`
`conda deactivate`

## Step 4 clone the artic-ncov2019 pipeline
`git clone --recursive https://github.com/artic-network/artic-ncov2019.git`
 
change into the repo directory
`cd artic-ncov2019`
]
## Step 5 create the conda environment for artic-ncov2019
`conda env create -n artic-ncov2019`
`conda activate artic-ncov2019`
`conda config --set channel_priority false`
`conda install artic-network::rampart=1.2.0`
`conda install snakemake-minimal=5.8.1`
`conda install -c bioconda -c conda-forge artic`

## Step 6 replace the primer scheme directory under artic-ncov2019 
`rm -r primer_schemes` 
`cd ..`

copy the primer_schemes_artic-ncov2019 folder from target folder to the artic-ncov2019 folder and rename as primer_schemes 
`cp -r ./target/primer_schemes_artic-ncov2019 ./artic-ncov2019/primer_schemes`

# Running the pipeline:

## activate the environment

`conda activate target`

## process the raw data to sample consensus sequences

`python target.py -in <path to fastq/fast5 folders> -s < path and name of sample_names.csv file>`
