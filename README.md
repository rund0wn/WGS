# MetaFlow
A Nextflow workflow to process short-read metagenome data. The workflow includes the following steps:
1. QC: fastp and fastQC
2. Assembly: MEGAHIT and metaQUAST
3. Annotation: prokka
4. Taxonomic assignment: metaphlan

#### Notes:
- The workflow is configured to work with docker or singularity. The singularity profile works with SLURM by default, an sbatch job can be submitted with the available example script.
- Fastq files must be named *_L001_R{1,2}_001.fastq.gz
![Screenshot from 2024-05-02 12-52-47](https://github.com/rund0wn/WGS/assets/107937921/76b46983-35c8-45f1-80f9-60d84c42b088)

## Setup:
- Install Nextflow
- Install Singularity or Docker
- Create 'data' folder in main directory and change relevant paths in the [config file](nextflow.config) (lines 13 and 24)

## To run:
nextflow run WGS.nf --in_dir directory/with/fastq/files -profile (docker OR singularity)

##### Workflow can be tested with:
nextflow run WGS.nf --in_dir test_reads -profile (docker OR singularity)
