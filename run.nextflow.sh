#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

j=${1:-'0'}
snapshot_date=${2:-'2022-02-25'}
batch_size=${3:-'1'}
pipeline=${4:-'nanopore'}
nextflow_script=${5:-"${HOME}/ENA_SARS_Cov2_nanopore/main.nf"}
profile=${6:-'gls'}
root_dir=${7:-'gs://prj-int-dev-covid19-nf-gls'}
dataset_name=${8:-'sarscov2_metadata'}
project_id=${9:-'prj-int-dev-covid19-nf-gls'}

#offset=$((j * batch_size))
table_name="${pipeline}_to_be_processed"
config_dir=$(dirname "${nextflow_script}")

# OFFSET ${offset}
sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size}"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
  | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv"

# Reserve ${table_name}_${j}.tsv
gsutil -m cp "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" "gs://${dataset_name}/${table_name}_${j}.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=0 --field_delimiter=tab \
  --autodetect --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${table_name}_${j}.tsv"

nextflow -C "${config_dir}/nextflow.config" run "${nextflow_script}" -profile "${profile}" \
      --INDEX "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" \
      --OUTDIR "${root_dir}/${snapshot_date}/${pipeline}_${j}/publishDir" \
      --STOREDIR "${root_dir}/${snapshot_date}/${pipeline}_${j}/storeDir" \
      -w "${root_dir}/${snapshot_date}/${pipeline}_${j}/workDir" \
      -with-tower

"${DIR}/update.metadata.sh" "${j}" "${snapshot_date}" "${pipeline}" "${dataset_name}" "${project_id}"
"${DIR}/update.receipt.sh" "${j}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
"${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"
