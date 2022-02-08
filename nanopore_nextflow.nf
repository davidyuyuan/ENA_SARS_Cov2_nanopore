#!/usr/bin/env nextflow

params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"

params.INDEX = "gs://prj-int-dev-covid19-nf-gls/prepro/nanopore.index.tsv"
params.SECRETS = "gs://prj-int-dev-covid19-nf-gls/prepro/projects_accounts.csv"

params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/storeDir"
params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/results"

//import nextflow.splitter.CsvSplitter
nextflow.enable.dsl=2

//def fetchRunAccessions(String tsv ) {
//    CsvSplitter splitter = new CsvSplitter().options( header:true, sep:'\t' )
//    BufferedReader reader = new BufferedReader( new FileReader( tsv ) )
//    splitter.parseHeader( reader )
//
//    List<String> run_accessions = [] as List<String>
//    Map<String,String> row
//    while( row = splitter.fetchRecord( reader ) ) {
//        run_accessions.add( row['run_accession'] )
//    }
//    return run_accessions
//}

process map_to_reference {
//    publishDir params.OUTDIR, mode:'copy'
    storeDir params.STOREDIR

    cpus 4 /* more is better, parallelizes very well*/
    memory '8 GB'
    container 'davidyuyuan/ena-sars-cov2-nanopore'  //:1.0'

    input:
    tuple(val(sampleId), file(input_file))
    path(sars2_fasta)
    path(sars2_fasta_fai)

    output:
//    file("${sampleId}_output/${sampleId}.bam")
//    file("${sampleId}_output/${sampleId}.coverage.gz")
//    file("${sampleId}_output/${sampleId}.annot.vcf.gz")
//    file("${sampleId}_output/${sampleId}_filtered.vcf.gz")
//    file("${sampleId}_output/${sampleId}_consensus.fasta.gz")
    val(sampleId, emit: sample_id)
    file("${sampleId}_output.tar.gz", emit: output_tgz)
    file("${sampleId}_filtered.vcf.gz", emit: filtered_vcf_gz)
    file("${sampleId}_consensus.fasta.gz", emit: consensus_fasta_gz)

    script:
    // curl -o ${sampleId}_1.fastq.gz \$(cat ${input_file})
    """
    wget -t 0 -O ${sampleId}_1.fastq.gz \$(cat ${input_file})
    cutadapt -u 30 -u -30 -o ${sampleId}.trimmed.fastq ${sampleId}_1.fastq.gz -m 75 -j ${task.cpus} --quiet

    minimap2 -Y -t ${task.cpus} -x map-ont -a ${sars2_fasta} ${sampleId}.trimmed.fastq | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${sampleId}.bam
    samtools index -@ ${task.cpus} ${sampleId}.bam

    bam_to_vcf.py -b ${sampleId}.bam -r ${sars2_fasta} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${sampleId}.vcf
    filtervcf.py -i ${sampleId}.vcf -o ${sampleId}_filtered.vcf
    bgzip ${sampleId}_filtered.vcf

    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${sampleId}.bam > ${sampleId}.pileup
    cat ${sampleId}.pileup | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
    tabix ${sampleId}_filtered.vcf.gz
    vcf2consensus.py -v ${sampleId}_filtered.vcf.gz -d ${sampleId}.coverage -r ${sars2_fasta} -o ${sampleId}_consensus.fasta -dp 30 -n ${sampleId}
    bgzip ${sampleId}.coverage
    bgzip ${sampleId}_consensus.fasta

    java -Xmx4g -jar /opt/conda/share/snpeff-5.0-1/snpEff.jar -q -no-downstream -no-upstream -noStats NC_045512.2 ${sampleId}.vcf > ${sampleId}.annot.vcf
    bgzip ${sampleId}.vcf
    bgzip ${sampleId}.annot.vcf

    mkdir -p ${sampleId}_output
    mv ${sampleId}.bam ${sampleId}.coverage.gz ${sampleId}.annot.vcf.gz ${sampleId}_output
    tar -zcvf ${sampleId}_output.tar.gz ${sampleId}_output
    mv ${sampleId}_filtered.vcf.gz ${sampleId}_consensus.fasta.gz ${sampleId}_output
    """
}

process ena-analysis-submit {
    publishDir params.OUTDIR, mode:'copy'
    storeDir params.STOREDIR
    secret 'webin_password'

    container 'davidyuyuan/ena-sars-cov2-nanopore:1.0'

    input:
    val(sampleId)
    file(output_tgz)
    file(filtered_vcf_gz)
    file(consensus_fasta_gz)
    path(projects_accounts_csv)

    output:
    file("${sampleId}_output.tar.gz")
    file("${sampleId}_output/${sampleId}_filtered.vcf.gz")
    file("${sampleId}_output/${sampleId}_consensus.fasta.gz")

    script:
    """
    cp -f config.yaml /usr/local/bin/config.yaml

    line=\$(grep "PRJEB43947" "${projects_accounts_csv}")
    webin_id=\$(echo "${line}" | cut -d ',' -f 4)
    webin_password=\$(echo "${line}" | cut -d ',' -f 5)

    analysis_submission.py -p PRJEB43947 -r ${sampleId} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au ${webin_id} -ap ${webin_password} -t
    analysis_submission.py -p PRJEB45554 -r ${sampleId} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au ${webin_id} -ap ${webin_password} -t
    analysis_submission.py -p PRJEB45619 -r ${sampleId} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS_VCF -au ${webin_id} -ap ${webin_password} -t

    mv ${output_tgz} ${sampleId}_output.tar.gz
    mv ${filtered_vcf_gz} ${sampleId}_filtered.vcf.gz
    mv ${consensus_fasta_gz} ${sampleId}_consensus.fasta.gz
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
    ena-analysis-submit(
            map_to_reference.out.sample_id,
            map_to_reference.out.output_tgz,
            map_to_reference.out.filtered_vcf_gz,
            map_to_reference.out.consensus_fasta_gz,
            params.SECRETS)
}
