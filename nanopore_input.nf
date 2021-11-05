#!/usr/bin/env nextflow

params.OUTDIR = "results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.index = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"

Channel
    .fromPath(params.index)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.run_accession, file('ftp://'+row.fastq_ftp)) }
    .set { samples_ch }

process check_input {
    input:
    set sampleId, input_file from samples_ch

    script:
    println("Count Fastq: " + input_file.countFastq())
    """
    echo "Run accession: $sampleId"
    echo "Input file: $input_file"
    """
}
