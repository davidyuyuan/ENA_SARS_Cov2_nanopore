#!/usr/bin/env nextflow

params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"
params.INDEX = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"
params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/storeDir"

import nextflow.splitter.CsvSplitter
nextflow.enable.dsl=2


def fetchRunAccessions(String tsv ) {
    CsvSplitter splitter = new CsvSplitter().options( header:true, sep:'\t' )
    BufferedReader reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )

    List<String> run_accessions = [] as List<String>
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
        run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

process download_fastq {
    storeDir params.STOREDIR

    // Use GLS default 1 CPU 1 GB and default quay.io/nextflow/bash
    // cpus 2
    // memory '1 GB'

    input:
    tuple val(sampleId), file(input_file)
    output:
    tuple val(sampleId), file("${sampleId}_1.fastq.gz")

    script:
    // curl -o ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
    wget -t 0 -O ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
}

process cut_adapters {
    storeDir params.STOREDIR

    // Use GLS default 1 CPU 1 GB
    // cpus 2
    // memory '1 GB'
    container 'kfdrc/cutadapt'

    input:
    tuple val(sampleId), file(input_file)

    output:
    tuple val(sampleId), file("${sampleId}.trimmed.fastq")

    script:
    """
    cutadapt -u 30 -u -30 -o ${sampleId}.trimmed.fastq ${input_file} -m 75 -j ${task.cpus} --quiet
    """
}

process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'
    storeDir params.STOREDIR

    cpus 8 /* more is better, parallelizes very well*/
    memory { 8.GB * task.attempt } //'8 GB'
    container 'davidyuyuan/ena-sars-cov2-nanopore'

    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(sampleId), file(trimmed)
    path(sars2_fasta)
    path(sars2_fasta_fai)

    output:
    file("${sampleId}.bam")
    file("${sampleId}_filtered.vcf.gz")
    file("${sampleId}.coverage")
    file("${sampleId}_consensus.fasta.gz")
    file("${sampleId}.annot.vcf")

    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${sars2_fasta} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${sampleId}.bam
    samtools index -@ ${task.cpus} ${sampleId}.bam
    bam_to_vcf.py -b ${sampleId}.bam -r ${sars2_fasta} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${sampleId}.vcf
    filtervcf.py -i ${sampleId}.vcf -o ${sampleId}_filtered.vcf
    bgzip ${sampleId}.vcf
    bgzip ${sampleId}_filtered.vcf
    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${sampleId}.bam > ${sampleId}.pileup
    cat ${sampleId}.pileup | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
    tabix ${sampleId}.vcf.gz
    vcf2consensus.py -v ${sampleId}.vcf.gz -d ${sampleId}.coverage -r ${sars2_fasta} -o ${sampleId}_consensus.fasta -dp 30 -n ${sampleId}
    bgzip ${sampleId}_consensus.fasta
    zcat ${sampleId}.vcf.gz | sed "s/^NC_045512.2/NC_045512/" > ${sampleId}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
}

workflow {
//    Requires local input.
//    accessions = fetchRunAccessions(params.INDEX)
//    data = Channel.fromSRA( accessions )
//    data.view()
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header:true, sep:'\t')
            .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp) }

    download_fastq(data)
    cut_adapters(download_fastq.out)
    map_to_reference(cut_adapters.out, params.SARS2_FA, params.SARS2_FA_FAI)
}
