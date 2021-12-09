#!/usr/bin/env bash

## DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

snapshot_date=${1:-'2021-12-13'}
pipeline=${2:-'nanopore'}
nextflow_script=${3:-"${HOME}/ENA_SARS_Cov2_nanopore/nanopore_nextflow.nf"}
concurrent_runs=${4:-'2'}
batch_size=${5:-'10000'}
dataset_name=${6:-'datahub_metadata'}
project_id=${7:-'prj-int-dev-covid19-nf-gls'}

# Result directories
output_dir="${DIR}/results/${snapshot_date}/output"
filteredvcf_dir="${DIR}/results/${snapshot_date}/filteredvcf"
consensus_dir="${DIR}/results/${snapshot_date}/consensus"
mkdir -p "${output_dir}" "${filteredvcf_dir}" "${consensus_dir}"

# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
batches=$(( row_count / batch_size + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}."

# Results and metadata
function gen_metadata {
  local input_file=$1
  local index_tsv=$2
  local to_be_submitted=$3
  local metadata=$4

  true > "${to_be_submitted}"
  true > "${metadata}"
  while IFS="" read -r f || [ -n "$f" ]
  do
    filename=$(basename "${f}")
    run_accession=$(echo "${filename}" | cut -d '_' -f 1)

    sample_accession=$(grep "${run_accession}" "${index_tsv}" | cut -f1-3)
    printf '%s\t%s\t%s\n' "${run_accession}" "${f}" "${sample_accession}" >> "${to_be_submitted}"

    sample_accession=$(grep "${run_accession}" "${index_tsv}" | cut -f3-4,10-13)
    printf '%s\t%s\n' "${run_accession}" "${sample_accession}" >> "${metadata}"
  done < "${input_file}"
}

for ((i=0; i<batches; i+=concurrent_runs)); do
  for ((j=i; j<i+concurrent_runs&&j<batches; j++)); do
    offset=$((j * batch_size))
    echo "Processing from offset ${offset}..."

    sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} ORDER BY run_accession DESC LIMIT ${batch_size} OFFSET ${offset}"
    bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
      | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv"

    echo "${HOME}/gcp-nf/gls/nextflow -C ${HOME}/gcp-nf/gls/nextflow.config run ${nextflow_script} -w gs://${project_id}/${snapshot_date}/${pipeline}_${j}/workDir --INDEX ${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv --OUTDIR gs://${project_id}/${snapshot_date}/${pipeline}_${j}/results --STOREDIR gs://${project_id}/${snapshot_date}/${pipeline}_${j}/storeDir -profile gls --resume -with-tower &"
    "${HOME}/gcp-nf/gls/nextflow" -C "${HOME}/gcp-nf/gls/nextflow.config" run "${nextflow_script}" \
      -w "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/workDir" \
      --INDEX "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" \
      --OUTDIR "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/results" \
      --STOREDIR "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/storeDir" \
      -profile gls --resume -with-tower &
  done
  wait

  for ((j=i; j<i+concurrent_runs&&j<batches; j++)); do
    echo "Lists of analysis URLs: ${j}..."
    gsutil -m ls "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/results/*_output.tar.gz" > "${output_dir}/${pipeline}_urls_${j}.txt" &
    gsutil -m ls "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/results/*/*_filtered.vcf.gz" > "${filteredvcf_dir}/${pipeline}_urls_${j}.txt" &
    gsutil -m ls "gs://${project_id}/${snapshot_date}/${pipeline}_${j}/results/*/*_consensus.fasta.gz" > "${consensus_dir}/${pipeline}_urls_${j}.txt" &
  done
  wait

  for ((j=i; j<i+concurrent_runs&&j<batches; j++)); do
    echo "Lists of metadata: ${j}..."
    gen_metadata "${output_dir}/${pipeline}_urls_${j}.txt" "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" \
      "${output_dir}/${pipeline}_to_be_submitted_${j}.tsv" "${output_dir}/${pipeline}_metadata_${j}.tsv" &
    gen_metadata "${filteredvcf_dir}/${pipeline}_urls_${j}.txt" "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" \
      "${filteredvcf_dir}/${pipeline}_to_be_submitted_${j}.tsv" /dev/null &
    gen_metadata "${consensus_dir}/${pipeline}_urls_${j}.txt" "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" \
      "${consensus_dir}/${pipeline}_to_be_submitted_${j}.tsv" /dev/null &
  done
  wait
done
