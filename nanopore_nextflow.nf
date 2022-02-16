#!/usr/bin/env nextflow

params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"

params.INDEX = "gs://prj-int-dev-covid19-nf-gls/prepro/nanopore.index.tsv"
params.SECRETS = "gs://prj-int-dev-covid19-nf-gls/prepro/projects_accounts.csv"
//params.CONFIG = "gs://prj-int-dev-covid19-nf-gls/prepro/config.nanopore.yaml"

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
    container 'davidyuyuan/ena-sars-cov2-nanopore:1.0'

    input:
//    tuple val(run_accession), val(sample_accession), file(input_file)
    tuple val(run_accession), file(input_file)
    path(sars2_fasta)
    path(sars2_fasta_fai)

    output:
//    file("${run_accession}_output/${run_accession}.bam")
//    file("${run_accession}_output/${run_accession}.coverage.gz")
//    file("${run_accession}_output/${run_accession}.annot.vcf.gz")
//    file("${run_accession}_output/${run_accession}_filtered.vcf.gz")
//    file("${run_accession}_output/${run_accession}_consensus.fasta.gz")
    //    , emit: sample_id//    , emit: output_tgz//    , emit: filtered_vcf_gz//    , emit: consensus_fasta_gz
    val(run_accession)
//    val(sample_accession)
    file("${run_accession}_output.tar.gz")
    file("${run_accession}_filtered.vcf.gz")
    file("${run_accession}_consensus.fasta.gz")

    script:
    // curl -o ${run_accession}_1.fastq.gz \$(cat ${input_file})
    """
    wget -t 0 -O ${run_accession}_1.fastq.gz \$(cat ${input_file})
    cutadapt -u 30 -u -30 -o ${run_accession}.trimmed.fastq ${run_accession}_1.fastq.gz -m 75 -j ${task.cpus} --quiet

    minimap2 -Y -t ${task.cpus} -x map-ont -a ${sars2_fasta} ${run_accession}.trimmed.fastq | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${run_accession}.bam
    samtools index -@ ${task.cpus} ${run_accession}.bam

    bam_to_vcf.py -b ${run_accession}.bam -r ${sars2_fasta} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${run_accession}.vcf
    filtervcf.py -i ${run_accession}.vcf -o ${run_accession}_filtered.vcf
    bgzip ${run_accession}_filtered.vcf

    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${run_accession}.bam > ${run_accession}.pileup
    cat ${run_accession}.pileup | awk '{print \$2,","\$3,","\$4}' > ${run_accession}.coverage
    tabix ${run_accession}_filtered.vcf.gz
    vcf2consensus.py -v ${run_accession}_filtered.vcf.gz -d ${run_accession}.coverage -r ${sars2_fasta} -o ${run_accession}_consensus.fasta -dp 30 -n ${run_accession}
    bgzip ${run_accession}.coverage
    bgzip ${run_accession}_consensus.fasta

    java -Xmx4g -jar /opt/conda/share/snpeff-5.0-1/snpEff.jar -q -no-downstream -no-upstream -noStats NC_045512.2 ${run_accession}.vcf > ${run_accession}.annot.vcf
    bgzip ${run_accession}.vcf
    bgzip ${run_accession}.annot.vcf

    mkdir -p ${run_accession}_output
    mv ${run_accession}.bam ${run_accession}.coverage.gz ${run_accession}.annot.vcf.gz ${run_accession}_output
    tar -zcvf ${run_accession}_output.tar.gz ${run_accession}_output
    """
}

process ena_analysis_submit {
    publishDir params.OUTDIR, mode:'copy'
    storeDir params.STOREDIR

    container 'davidyuyuan/ena-analysis-submitter:2.0'

    input:
    val(run_accession)
//    val(sample_accession)
    file(output_tgz)
    file(filtered_vcf_gz)
    file(consensus_fasta_gz)
    path(projects_accounts_csv)
//    path(config_yaml)

    output:
    file("PRJEB43947/${run_accession}_output.tar.gz")
    file("PRJEB45554/${run_accession}_filtered.vcf.gz")
    file("PRJEB45619/${run_accession}_consensus.fasta.gz")
    file("successful_submissions.txt")

    script:
//    cp -f ${config_yaml} /usr/local/bin/config.yaml
//    cat /usr/local/bin/config.yaml
//    echo \${PWD}
//    cd /usr/local/bin/ && ./analysis_submission.py -t -s ${sample_accession} -p PRJEB43947 -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password} &
//    cd /usr/local/bin/ && ./analysis_submission.py -t -s ${sample_accession} -p PRJEB45554 -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password} &
//    cd /usr/local/bin/ && ./analysis_submission.py -t -s ${sample_accession} -p PRJEB45619 -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password} &
//    wait && mv /usr/local/bin/successful_submissions.txt successful_submissions.txt
//    mkdir -p ${run_accession}_output
//    mv ${output_tgz} ${filtered_vcf_gz} ${consensus_fasta_gz} ${run_accession}_output
    """
    echo ${projectDir}
    echo $launchDir
    echo $workDir
    echo $homeDir

    webin_line="\$(grep PRJEB43947 ${projects_accounts_csv})"
    webin_id="\$(echo \${webin_line} | cut -d ',' -f 4)"
    webin_password="\$(echo \${webin_line} | cut -d ',' -f 5)"

    analysis_submission.py -t -o . -p PRJEB43947 -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password} &
    analysis_submission.py -t -o . -p PRJEB45554 -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password} &
    analysis_submission.py -t -o . -p PRJEB45619 -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password} &
    wait
    
    mkdir -p PRJEB43947 && mv ${output_tgz} PRJEB43947
    mkdir -p PRJEB45554 && mv ${filtered_vcf_gz} PRJEB45554
    mkdir -p PRJEB45619 && mv ${consensus_fasta_gz} PRJEB45619
    """
}

workflow {
//    Requires local input.
//    accessions = fetchRunAccessions(params.INDEX)
//    data = Channel.fromSRA( accessions )
//    data.view()
//    , row.sample_accession
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header:true, sep:'\t')
            .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp) }

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI)
    ena_analysis_submit(map_to_reference.out, params.SECRETS)
//    , params.CONFIG)
}
