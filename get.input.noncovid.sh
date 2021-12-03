#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d 'result=read_run&query=tax_tree(11118)%20AND%20((%20instrument_platform%3D%22ILLUMINA%22%20))&fields=accession%2Csample_accession%2Cstudy_accession%2Ccollection_date%2Ccollection_date_submitted%2Ccountry%2Cfirst_created%2Cinstrument_platform%2Cinstrument_model%2Ctax_id%2Cstrain%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_md5&excludeAccessionType=taxon&excludeAccessions=2697049&format=tsv' \
  "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/illumina.noncovid.index.tsv"

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d 'result=read_run&query=tax_tree(11118)%20AND%20((%20instrument_platform%3D%22OXFORD_NANOPORE%22%20))&fields=accession%2Csample_accession%2Cstudy_accession%2Ccollection_date%2Ccollection_date_submitted%2Ccountry%2Cfirst_created%2Cinstrument_platform%2Cinstrument_model%2Ctax_id%2Cstrain%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_md5&excludeAccessionType=taxon&excludeAccessions=2697049&format=tsv' \
  "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/nanopore.noncovid.index.tsv"

nohup ~/gcp-nf/gls/nextflow -C ~/gcp-nf/gls/nextflow.config run ~/ENA_SARS_Cov2_nanopore/nanopore_nextflow.nf \
  -w gs://prj-int-dev-covid19-nf-gls/nanopore/noncovid --INDEX ~/ENA_SARS_Cov2_nanopore/nanopore.noncovid.index.tsv \
  --STOREDIR "gs://prj-int-dev-covid19-nf-gls/noncovid/storeDir" --OUTDIR "gs://prj-int-dev-covid19-nf-gls/noncovid/results" \
  -profile gls --resume -with-tower
