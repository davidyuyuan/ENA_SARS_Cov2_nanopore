#!/usr/bin/env bash

snapshot_date=${1:-'2021-12-13'}
batches=${2:-'3'}
batch_size=${3:-'100000'}
table_name=${4:-'nanopore_to_be_processed'}
dataset_name=${5:-'datahub_metadata'}
project_id=${6:-'prj-int-dev-covid19-nf-gls'}

## DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

SQL_COUNT="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"

row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${SQL_COUNT}" | grep -v total)
echo "Row count: ${row_count}."
for ((i=0; i<batches; i++))
do
  offset=$((i * batch_size))
  if [[ ${offset} -lt ${row_count} ]]; then
    sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} ORDER BY run_accession DESC LIMIT ${batch_size} OFFSET ${offset}"
    bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" > "${DIR}/results/${snapshot_date}/${table_name}_${i}.csv"

    ~/gcp-nf/gls/nextflow -C ~/gcp-nf/gls/nextflow.config run ~/ENA_SARS_Cov2_nanopore/nanopore_nextflow.nf \
      -w "gs://prj-int-dev-covid19-nf-gls/${snapshot_date}/nanopore_${i}/workDir" --INDEX "${DIR}/results/${snapshot_date}/${table_name}_${i}.csv" \
      --OUTDIR "gs://prj-int-dev-covid19-nf-gls/${snapshot_date}/nanopore_${i}/results" --STOREDIR "gs://prj-int-dev-covid19-nf-gls/${snapshot_date}/nanopore_${i}/storeDir" \
      -profile gls --resume -with-tower &
  fi
done
