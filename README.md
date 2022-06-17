# ClosTyper   

ClosTyper is a pipeline tool designed to automate characterization and genotyping of selected _Clostrdia_ species using the whole genome sequencing data   
ClosTyper is written in Snakemake that allows reproducibility and scalability of the intergrated workflow. <br><br>
&#9995; **This tool is under active development** &#9995;   
You may think of use WGSBAC (https://gitlab.com/FLI_Bioinfo/WGSBAC) that combine different generalized modules for bacterial characterization based on WGS data  

## Requirements 
####  Software requirements
before start, you need to make sure that the follwoing software are installed in your system and that they are available in your $PATH 
1.  `any2fasta`: Convert various sequence formats to FASTA  (https://github.com/tseemann/any2fasta)   
    - Installation instructions: https://github.com/tseemann/any2fasta#github 
    - Required version: 0.4.2 ro later   
    
2.  `snakemake`: a workflow management system that creates reproducible and scalable data analyses. Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment (https://snakemake.github.io/)
    - Installation instructions: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html 
    - Required version: 6.15.5 or later    

3. `pigz`: A parallel implementation of gzip for modern multi-processor, multi-core machines (https://zlib.net/pigz/)
    - Installation instructions: https://zlib.net/pigz/
    - Required version: 2.3.4 or later    

4. `conda`: package, dependency and environment management system (https://docs.conda.io/en/latest/) 
    - Installation instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

####   Databases requirements
1. (Mini)Kraken2 database: taxonomic profiling of sequencing data (https://ccb.jhu.edu/software/kraken2/)
    - Download page: https://ccb.jhu.edu/software/kraken2/downloads.shtml

## Installation
#### Install clostyper from source 
These instructions will install the latest version of `ClosTyper`:

```
git clone https://gitlab.com/FLI_Bioinfo/clostyper.git  
ln -s `pwd`/clostyper/bin/clostyper /usr/local/bin # choose a folder in your $PATH

```  
This should be the directory structure of clostyper   

```
|-- bin
|   `-- clostyper
|-- config
|   `-- README.md
|-- dbs
|   |-- custom_dbs
|   `-- trimmomatic.fa
|-- README.md
`-- workflow
    |-- envs
    |-- rules
    |-- scripts
    `-- Snakefile
```

## Usage
### First execution 

On first execution of clostyper with `clostyper -h`, a config file will be automatically created <clostyper_config.txt> under the folder <clostyper/config>.   
This file should include the paths to databases and schemes to be used with `clostyper`.    
 &#9997; Important: at least the full path to the Kraken2 database must be provided.

## Basic usage [with test data]

#### download test data: 
```
wget -O test_data_clostyper.tar.gz https://zenodo.org/record/6656045/files/test_data_clostyper.tar.gz?download=1 
tar xzvf test_data_clostyper.tar.gz 
```
### Run clostyper 
* **To check raw data quality**  
This will only execute the fastp & kraken2    
Basic call: `clostyper --check_quality -d [FASTQ_DIRECTORY] -r [REFERENCE] [-o WORKING_DIRECTORY] `      
With the test data    
```
clostyper --check_quality -d test_data/input_data/ -r test_data/input_data/ref.gbk -o results_dir --kraken2db </path/2/kraken_database/> -A -T <cpus> 
```

* **To configure and execute the full pipeline**    
  &#9997; [WARNING: Many features of clostyper are still experimental. Please report issues to the issue tracker]  
This will only execute the full pipeline    
Basic call: `clostyper -d FASTQ/FASTA_DIRECTORY -r REFERENCE [-o WORKING_DIRECTORY] [-s SPECIES]`     
With the test data     
```
clostyper -d test_data/input_data/ -r test_data/input_data/ref.gbk -o results_dir --kraken2db </path/2/kraken_database/> -A -T <cpus> 
```
* **To search _C. difficile_ for mobile elements**   
To ONLY search your genomes for Clostridium difficile mobile elements, invoke the flag `--run_only_species_wf` together with flag `-s`    
Example with the test dataset
```
clostyper -d test_data/input_data/ -r test_data/input_data/ref.gbk -o results_dir  --kraken2db </path/2/kraken_database/> -s cdifficile --run_only_species_wf -A -T <cpus>
```

### Full usage options 
```
> clostyper -h
ClosTyper: Clostridia characterization and typing pipeline

Version: 0.1-beta (available at: https://gitlab.com/FLI_Bioinfo/ClosTyper)
USAGE:
 ---------------------
     clostyper --check_quality -d FASTQ_DIRECTORY [-o WORKING_DIRECTORY]
     clostyper -d FASTQ/FASTA_DIRECTORY -r REFERENCE [-o WORKING_DIRECTORY] [-s SPECIES] [--run_cgmlst]
     clostyper -t SAMPLE_TABLE -r REFERENCE [-o WORKING_DIRECTORY] [-s SPECIES] [--run_cgmlst]
 ---------------------
INPUT:
   -d, --fastx-directory          DIR, a directory where fastq reads or assembled genomes are present. Required unless -t flag was used
                                    Format: [ID]_{1,2}.fastq{.gz} [ID]_S*_R{1,2}_001.fastq{.gz} [ID]_R{1,2}.fastq{.gz} OR [ID].{fasta,fna,fa}
   -t, --sample_table             FILE (tab delimited), a four-columns based table. See an example in the documentation!
                                    Required unless the -d flag was used. If -t and -d flags were activated, -d will be ignored
   -r, --reference                FILE, Reference genome. Format: {ID}.{gbk,fasta,gff,embl} (required)
   -s, --species                  Run species-specific workflow (default: False; run only the general workflow)
                                    Currently supported Clostridia species are: cdifficile
OUTPUT:
   -o, --output-directory         DIR, output directory for the snakemake results (default: output_dir_[timestamp]/)
   -w, --overwrite                Overwrite an existing directory with the results. Useful to append results to previous runs
   -q, --quiet                    Suppress clostyper messages. Report only warnings, errors and the snakemake call
WORKFLOW:
   -Q, --check_quality            Quickly perform quality assurance on Illumina data [recommended before doing analysis] [EXPERIMENTAL]
   --run_cgmlst                   Do cgMLST analyis using the chewiesnake pipeline [EXPERIMENTAL]
   --run_only_species_wf          Execute ONLY the specified species-specific workflow. Require '-s' [EXPERIMENTAL]
   --select_reference             Select appropriate reference for the SNP anlysis of the dataset [EXPERIMENTAL]
   --snp_pipeline                 Select which pipeline to call SNPs [EXPERIMENTAL]
                                    Supported SNP pipelines are: snippy, reddog, nasp, cfsanpipeline
   --kraken2db                    Path to (Mini)kraken2 DB
   --disable_pangenome            Disable the pangenome analyis [EXPERIMENTAL]
   --disable_report               Do not make the html report [EXPERIMENTAL]
OTHERS:
   -A, --autorun                  Automatically run snakemake workflow after configuration (default: False)
   -T, --threads                  Number of threads to use (default: 16)
   --check_dep                    Check if dependencies are ok and then exit
   --no-color                     Do not use a colored output (default: False)
HELP:
   -h, --help                     Show this help and exit
   --help_all                     Show extended help for all software settings options [EXPERIMENTAL]
   --version                      Show clostyper's version number and exit
   --citation                     Show clostyper's citation and exit

```
