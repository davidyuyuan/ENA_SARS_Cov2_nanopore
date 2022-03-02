#!/usr/bin/env bash

bad_records=${1:-'0'}
pipeline=${2:-'nanopore'}
dataset_name=${3:-'datahub_metadata'}
project_id=${4:-'prj-int-dev-covid19-nf-gls'}
location=${5:-'europe-west4'}

# Create bucket and dataset
gsutil ls -p "${project_id}" "gs://${dataset_name}" || gsutil mb -p "${project_id}" -l "${location}" "gs://${dataset_name}"
bq --project_id="${project_id}" show "${dataset_name}" || bq --location="${location}" mk --dataset "${project_id}:${dataset_name}"

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

############################################################
# Get all raw reads for Nanopore and Illumina into Big Query
############################################################
if [ "${pipeline}" = "nanopore" ]; then
  curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
    -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
    "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/${pipeline}_index.tsv"
fi
if [ "${pipeline}" = "illumina" ]; then
  curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
    -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22ILLUMINA%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
    "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/${pipeline}_index.tsv"
fi
gsutil -m cp "${DIR}/${pipeline}_index.tsv" "gs://${dataset_name}/${pipeline}_index.tsv"

# It is possible that ILLUMINA raw reads are mislabeled as OXFORD_NANOPORE: grep '_2.fastq.gz' nanopore_index.tsv | wc -l is 296
bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect --max_bad_records="${bad_records}" "${dataset_name}.${pipeline}_index" "gs://${dataset_name}/${pipeline}_index.tsv"

# copy table schema if it does not exist
sql="SELECT COUNT(*) AS total FROM ${dataset_name}.__TABLES_SUMMARY__ WHERE table_id = '""${pipeline}_processing'"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
if [ "${row_count}" = "0" ]; then
  sql="SELECT * FROM ${dataset_name}.${pipeline}_index LIMIT 0"
  bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --destination_table "${dataset_name}.${pipeline}_processing" "${sql}"
fi
