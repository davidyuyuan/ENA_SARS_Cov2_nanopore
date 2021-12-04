#!/usr/bin/env bash

staging_dir=${1:-'/mnt/result/from_gcs/nanopore/staging'}
date_submitted=${2:-'2021-11-29'}
#staging_dir='/Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/staging'
#date_submitted='2021-11-29'

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

index_tsv="${DIR}/../${date_submitted}/nanopore.index.tsv"
to_be_submitted="${DIR}/../${date_submitted}/output/to_be_submitted.tsv"

true > "${to_be_submitted}" || touch "${to_be_submitted}"
for f in "${staging_dir}"/*
do
  filename=$(basename "${f}")
  run_accession=$(echo "${filename}" | cut -d '_' -f 1)
  sample_accession=$(grep "${run_accession}" "${index_tsv}" | cut -f1-3)
  printf '%s\t%s\t%s\n' "${run_accession}" "${f}" "${sample_accession}" >> "${to_be_submitted}"
done
