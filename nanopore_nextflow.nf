#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.OUTDIR = "results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.index = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"

Channel
    .fromPath(params.index)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp) }
    .set { samples_ch }

process download_fastq {
    container 'davidyuyuan/samtools:dlf'
    // Use GLS default 1 CPU 1 GB
    // cpus 2
    // memory '1 GB'

    input:
    tuple sampleId, file(input_file) from samples_ch
    output:
    tuple sampleId, file("${sampleId}_1.fastq.gz") into fastq_ch

    script:
    // wget -O ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
    curl -o ${sampleId}_1.fastq.gz \$(cat ${input_file})
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

/*
 * Map reads to SARS-CoV-2 reference genome
 */ 
/*
process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'

    cpus 10 // more is better, parallelizes very well
    memory '10 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path trimmed from trimmed_ch
    path ref from params.SARS2_FA
    val run_id from params.RUN_ID
    
    output:
    path "${run_id}.bam" into sars2_aligned_reads_ch, sars2_aligned_reads_ch2
    path("${run_id}.bam")
    
    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${ref} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${run_id}.bam
    """
}

process check_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '4 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from sars2_aligned_reads_ch
    val run_id from params.RUN_ID
    path sars2_fasta from params.SARS2_FA

    output:
    path "${run_id}.pileup" into check_coverage_ch
    path("${run_id}.pileup")

    script:
    """
    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${bam} > \
    ${run_id}.pileup
    """
}

process make_small_file_with_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '4 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path pileup from check_coverage_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.coverage")

    script:
    """
    cat ${pileup} | awk '{print \$2,","\$3,","\$4}' > ${run_id}.coverage
    """
}

process bam_to_vcf {
    tag '$run_id'
    publishDir params.OUTDIR, mode:'copy'
    cpus 10
    memory '10 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'

    input:
    path bam from sars2_aligned_reads_ch2
    path ref from params.SARS2_FA
    val run_id from params.RUN_ID
    
    output:
    path "${run_id}.vcf" into vcf_ch

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    bam_to_vcf.py -b ${bam} -r ${ref} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${run_id}.vcf
    """
}
*/
process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'

    cpus 8 /* more is better, parallelizes very well*/
    memory '8 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    tuple sampleId, file(trimmed) from trimmed_ch
    path(ref) from params.SARS2_FA
    // val run_id from params.RUN_ID
    
    output:
//    tuple sampleId, file("${sampleId}.bam") into sars2_aligned_reads_ch
    file("${sampleId}.bam")
    tuple sampleId, file("${sampleId}.vcf") into vcf_ch
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
/*
process check_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 4
    memory '4 GB'
    container 'davidyuyuan/samtools:dlf'    // 'alexeyebi/bowtie2_samtools' // staphb/samtools

    input:
    tuple sampleId, file(bam) from sars2_aligned_reads_ch
    // val run_id from params.RUN_ID
    path(sars2_fasta) from params.SARS2_FA

    output:
    file("${sampleId}.pileup")
    file("${sampleId}.coverage")

    script:
    """
    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${bam} > ${sampleId}.pileup
    cat ${sampleId}.pileup | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
    """
}
*/
process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 8
    memory '8 GB'
    container 'alexeyebi/snpeff'    // nfcore/snpeff

    input:
    tuple sampleId, file(vcf) from vcf_ch
    // val run_id from params.RUN_ID

    output:
    file("${sampleId}.annot.vcf")

    script:
//    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
    cat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > ${sampleId}.newchr.vcf
    java -Xmx4g -jar /home/biodocker/bin/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
}
