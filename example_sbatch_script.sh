#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J WGS_NF
#SBATCH -o WGS_NF.%J.out
#SBATCH -e WGS_NF.%J.err
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=12

module load nextflow
module load singularity

nextflow run WGS.nf --in_dir directory/with/fastq/files -profile singularity
