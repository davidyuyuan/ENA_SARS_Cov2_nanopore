#!/usr/bin/env nextflow

//21.07.0-edge
//NXF_VER='21.10.1' nextflow run 'https://github.com/davidyuyuan/ENA_SARS_Cov2_nanopore/blob/master/nanopore_nextflow_ncbi.nf'

params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"

//params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/results"
//params.INDEX = "gs://prj-int-dev-covid19-nf-gls/data/nanopore.index.tsv"
//params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/storeDir"
params.OUTDIR = "gs://ncbi-sra-covid-input/results"
//params.INDEX = "gs://ncbi-sra-covid-input/data/nanopore.index.short.tsv"
params.STOREDIR = "gs://ncbi-sra-covid-input/storeDir"

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

workflow {
    accessions = fetchRunAccessions(params.INDEX)
    Channel.fromSRA(accessions, protocol:'https').view()
}
