#!/usr/bin/env bash

# Umbrella study for all COVID-19 submissions: https://www.ebi.ac.uk/ena/browser/view/PRJEB45555?show=component-projects
# Specific project for output.tar.gz: https://www.ebi.ac.uk/ena/browser/view/PRJEB43947

home_dir=${1:-'/mnt/result/from_gcs/nanopore/new'}
cd "${home_dir}" || exit
#home_dir='/Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/new'
#cd "${home_dir}" || exit
#scp -i ~/.ssh/david-yuan-tsi.pem centos@45.88.81.197:/mnt/result/from_gcs/nanopore/new/ERR4080473* ${home_dir}
#scp -i ~/.ssh/david-yuan-tsi.pem centos@45.88.81.197:/mnt/result/from_gcs/nanopore/new/ERR4080474* ${home_dir}
#scp -i ~/.ssh/david-yuan-tsi.pem centos@45.88.81.197:/mnt/result/from_gcs/nanopore/new/ERR4080475* ${home_dir}

for f in "${home_dir}"/*
do
  filename=$(basename "${f}")
  run_accession=$(echo "${filename}" | cut -d '.' -f 1 | cut -d '_' -f 1)
  mkdir -p "${home_dir}/${run_accession}_output"

  case "${filename}" in
  *.vcf | *.coverage )
    gzip "${home_dir}/${filename}"
    mv "${home_dir}/${filename}.gz" "${home_dir}/${run_accession}_output"
    ;;
  * )
    mv "${home_dir}/${filename}" "${home_dir}/${run_accession}_output"
    ;;
  esac
done

staging_dir="${home_dir}/../staging" && mkdir -p ${staging_dir}
for f in "${home_dir}"/*
do
  filename=$(basename "${f}")
  run_accession=$(echo "${filename}" | cut -d '_' -f 1)
  output_tgz=${run_accession}_output.tar.gz
  tar -zcvf "${home_dir}/${output_tgz}" "${run_accession}_output" && rm -R "${run_accession}_output" && mv "${home_dir}/${output_tgz}" "${staging_dir}"
#  printf '%s\n' "${run_accession} ${filename}"
done
