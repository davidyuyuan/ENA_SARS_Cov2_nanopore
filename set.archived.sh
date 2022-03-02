#!/usr/bin/env bash

snapshot_date=${1:-'2022-02-25'}
dataset_name=${2:-'sarscov2_metadata'}
project_id=${3:-'prj-int-dev-covid19-nf-gls'}

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

output_dir="${DIR}/results/${snapshot_date}/output"
mkdir -p "${output_dir}"

###########################################################
# Get all analyses archived under parent_study = PRJEB45555
###########################################################
## PRJEB43947 output_tgz
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=study_accession%3D%22PRJEB43947%22&fields=sample_accession%2Crun_ref&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > ./archived_output.tsv
## PRJEB45554 filtered_vcf_gz
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=study_accession%3D%22PRJEB45554%22&fields=sample_accession%2Crun_ref&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > ./archived_filtered_vcf.tsv
## PRJEB45619 consensus_fasta_gz
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=study_accession%3D%22PRJEB45619%22&fields=sample_accession%2Crun_ref&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > ./archived_consensus_fasta.tsv

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=parent_study%3D%22PRJEB45555%22&fields=sample_accession%2Crun_ref&format=tsv&limit=0' \
  "https://www.ebi.ac.uk/ena/portal/api/search" > "${output_dir}/analysis_archived.tsv"
gsutil -m cp "${output_dir}/analysis_archived.tsv" "gs://${dataset_name}/analysis_archived.tsv"
bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.analysis_archived" "gs://${dataset_name}/analysis_archived.tsv" \
  "analysis_accession:STRING,sample_accession:STRING,run_ref:STRING"
