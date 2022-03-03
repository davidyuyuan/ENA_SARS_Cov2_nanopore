#!/usr/bin/env nextflow

params.SARS2_FA = "/homes/sands/ref-seq/NC_045512.2.fa"
params.SARS2_FA_FAI = "/homes/sands/ref-seq/NC_045512.2.fa.fai"

params.INDEX = "/homes/sands/1k/nanopore_to_be_processed_1k.tsv"
params.STOREDIR = "/hps/nobackup/cochrane/ena/users/sands/1k/storeDir"
params.OUTDIR = "/hps/nobackup/cochrane/ena/users/sands/1k/results"

nextflow.enable.dsl=2

process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'
    storeDir params.STOREDIR

    cpus 4 /* more is better, parallelizes very well*/
    memory '8 GB'
    container 'davidyuyuan/ena-sars-cov2-nanopore:1.0'

    input:
    tuple val(sampleId), file(input_file)
    path(sars2_fasta)
    path(sars2_fasta_fai)

    output:
    file("${sampleId}_output.tar.gz")
    file("${sampleId}_output/${sampleId}_filtered.vcf.gz")
    file("${sampleId}_output/${sampleId}_consensus.fasta.gz")

    script:
    // curl -o ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
    wget -t 0 -O ${sampleId}_1.fastq.gz \$(cat ${input_file})
    cutadapt -u 30 -u -30 -o ${sampleId}.trimmed.fastq ${sampleId}_1.fastq.gz -m 75 -j ${task.cpus} --quiet

    minimap2 -Y -t ${task.cpus} -x map-ont -a ${sars2_fasta} ${sampleId}.trimmed.fastq | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${sampleId}.bam
    samtools index -@ ${task.cpus} ${sampleId}.bam

    bam_to_vcf.py -b ${sampleId}.bam -r ${sars2_fasta} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${sampleId}.vcf
    filtervcf.py -i ${sampleId}.vcf -o ${sampleId}_filtered.vcf
    bgzip ${sampleId}.vcf
    bgzip ${sampleId}_filtered.vcf

    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${sampleId}.bam > ${sampleId}.pileup
    cat ${sampleId}.pileup | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
    tabix ${sampleId}_filtered.vcf.gz
    vcf2consensus.py -v ${sampleId}_filtered.vcf.gz -d ${sampleId}.coverage -r ${sars2_fasta} -o ${sampleId}_consensus.fasta -dp 30 -n ${sampleId}
    bgzip ${sampleId}.coverage
    bgzip ${sampleId}_consensus.fasta

    zcat ${sampleId}.vcf.gz | sed "s/^NC_045512.2/NC_045512/" > ${sampleId}.newchr.vcf
    java -Xmx4g -jar /opt/conda/share/snpeff-5.0-1/snpEff.jar -q -no-downstream -no-upstream -noStats NC_045512.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    bgzip ${sampleId}.annot.vcf

    mkdir -p ${sampleId}_output
    mv ${sampleId}.bam ${sampleId}.coverage.gz ${sampleId}.annot.vcf.gz ${sampleId}_output
    tar -zcvf ${sampleId}_output.tar.gz ${sampleId}_output
    mv ${sampleId}_filtered.vcf.gz ${sampleId}_consensus.fasta.gz ${sampleId}_output
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

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI)
}