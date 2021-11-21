#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
  "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/nanopore.index.tsv"

NXF_VER=21.10.1 ~/gcp-nf/gls/nextflow -C ~/gcp-nf/gls/nextflow.config run ~/ENA_SARS_Cov2_nanopore/nanopore_nextflow.nf \
  -w gs://prj-int-dev-covid19-nf-gls/nanopore --INDEX ~/ENA_SARS_Cov2_nanopore/nanopore.index.short.tsv \
  -profile gls --resume -with-tower
