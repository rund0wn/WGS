#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rund0wn/WGS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/rund0wn/WGS
----------------------------------------------------------------------------------------
*/

params.in_dir = './' // Default input directory
params.out_dir = 'results' // Default output directory
params.skip_trimming = false
params.db_path = "/data/metaphlan_dirs"
// params.python_paths = "/data/"

/*
    ================================================================================
                                Preprocessing and QC
    ================================================================================
*/


process fastp {
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/trimmed_reads", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    def (read1, read2) = reads
    """
    fastp \
        -o ${pair_id}_L001_R1_001_trimmed.fastq.gz \
        -O ${pair_id}_L001_R2_001_trimmed.fastq.gz \
        -i ${read1} \
        -I ${read2}
    """
}

process FastQC {
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/fastQC", mode: 'copy'

    input:
    tuple val(pair_id), path(read_files)

    output:
    path "*.html", emit: fastqc_html
    path "*.zip", emit: fastqc_zip

    script:
    def (read1, read2) = read_files
    """
    fastqc ${read1} ${read2}
    """
}

/*
    ================================================================================
                                Taxonomic Assignment
    ================================================================================
*/

process AssignTaxa {
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/taxa", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)

    output:
    path "profiled_metagenome.txt", emit: assigned_taxa
    path "metagenome.bowtie2.bz2", emit: bowtie_taxa

    script:
    def (read1, read2) = trimmed_reads
    """
    metaphlan ${read1},${read2} --bowtie2out metagenome.bowtie2.bz2 --bowtie2db ${params.db_path} --nproc 5 --input_type fastq > profiled_metagenome.txt
    """
}

// // /*
// //     ================================================================================
// //                                      Assembly
// //     ================================================================================
// // */

process Assemble{
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/genome", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)

    output:
    tuple val(pair_id), file("${pair_id}_assembly/final.contigs.fa"), emit: assembly
    tuple val(pair_id), file("${pair_id}_assembly/checkpoints.txt")
    tuple val(pair_id), file("${pair_id}_assembly/done")
    tuple val(pair_id), file("${pair_id}_assembly/log")
    tuple val(pair_id), file("${pair_id}_assembly/intermediate_contigs")

    script:
    def (read1, read2) = trimmed_reads
    """
    megahit -1 ${read1} -2 ${read2} -o ${pair_id}_assembly
    """
}

process Evaluate{
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/genome/${pair_id}_assembly", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    path "quast_${pair_id}"

    script:
    """
    metaquast.py ${assembly} -o quast_${pair_id}
    """
}

process Bin{
    label 'metabat'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/genome/${pair_id}_assembly/bin", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    path "${pair_id}_bins/*.fa"

    script:
    """
    gzip -f ${assembly}
    metabat -i ${assembly}.gz -o ${pair_id}_bins/bin
    """
}

// // /*
// //     ================================================================================
// //                                      Annotation
// //     ================================================================================
// // */

// process PredictProteins{
//     label 'def'
//     publishDir "${params.out_dir}_results/genome/annotation", mode: 'copy'

//     input:
//     path file

//     output:
//     path "${file.baseName}.faa"

//     script:
//     """
//     prodigal -i ${file} -a ${file.baseName}.faa -p meta
//     """
// }

process AnnotateGenome {
    label 'def'
    tag "${pair_id}"
    publishDir "${params.out_dir}_results/genome/annotation", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    path "${pair_id}"

    script:
    """
    prokka --outdir ${pair_id} --prefix ${pair_id} ${assembly}
    """
}

workflow {
    read_pairs = Channel.fromFilePairs("${params.in_dir}/*_L001_R{1,2}_001.fastq.gz", size: 2)

    //Merge
    if (!params.skip_trimming) {
    fastp(read_pairs)
    FastQC(fastp.out.trimmed_reads)
    } else {
        FastQC(read_pairs)

    }

    Assemble(fastp.out.trimmed_reads)
    Evaluate(Assemble.out.assembly)
    // bins = Bin(Assemble.out.assembly).flatten()
    // AnnotateGenome(bins)
    AnnotateGenome(Assemble.out.assembly)
    AssignTaxa(fastp.out.trimmed_reads)
}

}
