#!/usr/bin/env nextflow

params.OUTDIR = "results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.index = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"

Channel
    .fromPath(params.index)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp) }
    .set { samples_ch }

process download_fastq {
    input:
    set sampleId, file(input_file) from samples_ch
    output:
    file "${sampleId}_1.fastq.gz" into fastq_ch

    script:
    """
    # /usr/bin/env groovy
    println("Count Fastq: " + input_file.countFastq())
    println("" + input_file.name)
    """
}



process check_input {
    input:
    set sampleId, file(input_file) from samples_ch
    output:
    file "${sampleId}_1.fastq.gz" into fastq_ch

    script:
 //   println("Count Fastq: " + input_file.countFastq())
 //   tmp_link = input_file.getName()
 //   input_file.copyTo(tmp_link)

 //   echo "Input file: $tmp_link"
    """
    echo "Run accession: $sampleId"
    cp $input_file ${sampleId}_1.fastq.gz    
    """
}

/*
process get_input {
    input:
    set sampleId, file(input_file) from samples_ch    
    output:
    file 'fname' into input_ch
//    Channel.fromPath(input_file).set{data_ch}

    script:    
//    fname = sampleId + '_1.fastq.gz'
//    input_file.copyTo(fname)
    """
    fname=${sampleId}_1.fastq.gz
    cp $input_file $fname
    """
}
*/
/*
Channel
    .from(samples_ch)
    .set {sampleId, input_file}
Channel
    .from()
params.path = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR196/002/SRR1965912/SRR1965912_1.fastq.gz"

Channel.fromPath(params.path).set{ch_data}

process test  {

  input:
  file f from ch_data

  script:
  """
  echo $f
  """

}
*/
