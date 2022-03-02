#!/usr/bin/env bash

dataset_name=${1:-'sarscov2_metadata'}
project_id=${2:-'prj-int-dev-covid19-nf-gls'}
location=${3:-'europe-west4'}

# Create bucket and dataset
gsutil ls -p "${project_id}" "gs://${dataset_name}" || gsutil mb -p "${project_id}" -l "${location}" "gs://${dataset_name}"
bq --project_id="${project_id}" show "${dataset_name}" || bq --location="${location}" mk --dataset "${project_id}:${dataset_name}"

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

############################################################
# Get all raw reads for Nanopore and Illumina into Big Query
############################################################
# It is possible that ILLUMINA raw reads are mislabeled as OXFORD_NANOPORE: grep '_2.fastq.gz' nanopore_index.tsv | wc -l is 296
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
#  -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22ILLUMINA%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/illumina_index.tsv" && \
#  gsutil -m cp "${DIR}/illumina_index.tsv" "gs://${dataset_name}/illumina_index.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect --max_bad_records=5 "${dataset_name}.sra_index" "gs://${dataset_name}/illumina_index.tsv"

#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
#  -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > "${DIR}/nanopore_index.tsv" && \
#  gsutil -m cp "${DIR}/nanopore_index.tsv" "gs://${dataset_name}/nanopore_index.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
  --max_bad_records=296 "${dataset_name}.sra_index" "gs://${dataset_name}/nanopore_index.tsv"

# copy table schema if it does not exist
sql="SELECT COUNT(*) AS total FROM ${dataset_name}.__TABLES_SUMMARY__ WHERE table_id = '""sra_processing'"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
if [ "${row_count}" = "0" ]; then
  sql="SELECT * FROM ${dataset_name}.sra_index LIMIT 0"
  bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --destination_table "${dataset_name}.sra_processing" "${sql}"
fi
