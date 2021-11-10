#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.index = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"

Channel
    .fromPath(params.index)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp) }
    .set { samples_ch }

process download_fastq {
    // Use GLS default 1 CPU 1 GB and default quay.io/nextflow/bash
    // cpus 2
    // memory '1 GB'

    input:
    tuple sampleId, file(input_file) from samples_ch
    output:
    tuple sampleId, file("${sampleId}_1.fastq.gz") into fastq_ch

    script:
    // curl -o ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
    wget -t 3 -O ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
}

/*
 * Trim 30 nucleotides of each end of the reads using cutadapt to ensure that primer derived sequences are not used to generate a consensus sequence
 */  
process cut_adapters {
    // Use GLS default 1 CPU 1 GB
    // cpus 2
    // memory '1 GB'
    container 'kfdrc/cutadapt'
    
    input:
    tuple sampleId, file(input_file) from fastq_ch
    
    output:
    tuple sampleId, file('trimmed.fastq') into trimmed_ch
    
    script:
    """
    cutadapt -u 30 -u -30 -o trimmed.fastq ${input_file} -m 75 -j ${task.cpus} --quiet
    """
}

process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'

    cpus 8 /* more is better, parallelizes very well*/
    memory '8 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    tuple sampleId, file(trimmed) from trimmed_ch
    path(ref) from params.SARS2_FA
    
    output:
    file("${sampleId}.bam")
    val sampleId into vcf_ch1
    path "${sampleId}.vcf" into vcf_ch2
    file("${sampleId}.pileup")
    file("${sampleId}.coverage")
    
    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${ref} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${sampleId}.bam
    samtools index -@ ${task.cpus} ${sampleId}.bam
    bam_to_vcf.py -b ${sampleId}.bam -r ${ref} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${sampleId}.vcf
    samtools mpileup -a -A -Q 0 -d 8000 -f ${ref} ${sampleId}.bam > ${sampleId}.pileup
    cat ${sampleId}.pileup | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
    """
}

process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 8
    memory '8 GB'
    container 'alexeyebi/snpeff'

    input:
    val sampleId from vcf_ch1
    path vcf from vcf_ch2

    output:
    file("${sampleId}.annot.vcf")
    file("${sampleId}.newchr.vcf")

    script:
//    java -Xmx4g -jar /home/biodocker/bin/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
    cat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > ${sampleId}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
}
