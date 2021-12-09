#!/usr/bin/env bash

# Umbrella study for all COVID-19 submissions: https://www.ebi.ac.uk/ena/browser/view/PRJEB45555?show=component-projects
# Specific project for output.tar.gz: https://www.ebi.ac.uk/ena/browser/view/PRJEB43947

## DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

snapshot_date=${1:-'2021-12-13'}
dataset_name=${2:-'datahub_metadata'}
project_id=${3:-'prj-int-dev-covid19-nf-gls'}

output_dir="${DIR}/results/${snapshot_date}/output/"
mkdir -p "${output_dir}"

#for ((i=0; i<batches; i++))
#do
#  offset=$((i * batch_size))
#  if [[ ${offset} -lt ${row_count} ]]; then
#    gsutil -m ls "gs://${project_id}/${snapshot_date}/nanopore_${i}/results/*" > "${output_dir}/nanopore_metadata.tsv"
#  fi
#done


#gsutil -m ls gs://${project_id}/${snapshot_date}/nanopore_${i}/results/* > all_results.txt

#gsutil -m cp "${output_dir}/new_nanopore_metadata.tsv" "gs://${dataset_name}/nanopore_metadata.tsv"
#
#bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
#  --autodetect "${dataset_name}.nanopore_metadata" "gs://${dataset_name}/nanopore_metadata.tsv"
