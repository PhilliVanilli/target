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

## Step 3.1 activate the conda env
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

python target/target.py --help  
usage: target.py [-h] -in PROJECT_PATH -r  
                 {ChikAsian_V1_400,ChikECSA_V1_800,ZikaAsian_V1_400,SARS2_V1_800,SARS2_V1_400,RSVA_V1_3000,RSVB_V1_3000,DENV1_V1_400,DENV2_V1_400}  
                 [-rs REFERENCE_START] [-re REFERENCE_END] [-mi MIN_LEN]  
                 [-ma MAX_LEN] [-d MIN_DEPTH] [--run_step RUN_STEP]
                 [--run_step_only] [-b {0,1}] [-m] [-a]
                 [-c {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}] [-g GPU_THREADS]  
                 [-gb GPU_BUFFERS] [--use_gaps] [--use_bwa] -p GUPPY_PATH  
                 [-rt]  

Process raw nanopore reads to fasta consensus sequences

optional arguments:
  -h, --help            show this help message and exit
  -in PROJECT_PATH, --project_path PROJECT_PATH
                        The path to the directory containing the 'fast5' and 'fastq' folders
  -r {ChikAsian_V1_400,ChikECSA_V1_800,ZikaAsian_V1_400,SARS2_V1_800,SARS2_V1_400,RSVA_V1_3000,RSVB_V1_3000,DENV1_V1_400,DENV2_V1_400}, --reference {ChikAsian_V1_400,ChikECSA_V1_800,ZikaAsian_V1_400,SARS2_V1_800,SARS2_V1_400,RSVA_V1_3000,RSVB_V1_3000,DENV1_V1_400,DENV2_V1_400}
                        The reference genome and primer scheme to use (default: None)
  -rs REFERENCE_START, --reference_start REFERENCE_START
                        The start coordinate of the reference sequence for read mapping (default: 1)
  -re REFERENCE_END, --reference_end REFERENCE_END
                        The end coordinate of the reference sequence for read mapping. Default = full length (default: False)
  -mi MIN_LEN, --min_len MIN_LEN
                        The minimum read length allowed:
                         = 300 for 400bp amplicon design
                         = 700 for 800bp amplicon design (default: 700)
  -ma MAX_LEN, --max_len MAX_LEN
                        The maximum read length allowed:
                         = 500 for 400bp amplicon design
                         = 900 for 800bp amplicon design (default: 900)
  -d MIN_DEPTH, --min_depth MIN_DEPTH
                        The minimum coverage to call a position in the MSA to consensus (default: 100)
  --run_step RUN_STEP   Run the pipeline starting at this step:
                        --run_step 0 = basecall reads with Guppy
                        --run_step 1 = demultiplex with Guppy
                        --run_step 2 = size filer and rename demultiplexed fastq file
                        --run_step 3 = concatenate demultiplexed files into sample files
                        --run_step 4 = run read mapping and all the variant calling steps on each sample
                         (default: 0)
  --run_step_only       Only run the step specified in 'run_step' (default: False)
  -b {0,1}, --basecall_mode {0,1}
                        0 = basecall in r10.4.1 kit14 mode
                        1 = basecall in r9.4.1 mode
                         (default: 1)
  -m, --msa             Generate consensus from MSA (default: False)
  -a, --art             Generate consensus with Artic pipeline (default: False)
  -c {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, --cpu_threads {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
                        The number of cpu threads to use for bwa, nanopolish etc... (default: 16)
  -g GPU_THREADS, --gpu_threads GPU_THREADS
                        The number of gpu threads to use ... (default: 8)
  -gb GPU_BUFFERS, --gpu_buffers GPU_BUFFERS
                        The number of gpu buffers to use for demultiplexing (default: 15)
  --use_gaps            use gap characters when making the consensus sequences (default: )
  --use_bwa             use bwa instead of minimap2 to map reads to reference (default: )
  -p GUPPY_PATH, --guppy_path GUPPY_PATH
                        The path to the guppy executables eg: '.../ont-guppy/bin/'
  -rt, --real_time      start basecalling fast5 files in batches during sequencing (default: False)
