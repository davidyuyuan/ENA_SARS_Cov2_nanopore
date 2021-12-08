#!/usr/bin/env bash

# Umbrella study for all COVID-19 submissions: https://www.ebi.ac.uk/ena/browser/view/PRJEB45555?show=component-projects
# Specific project for output.tar.gz: https://www.ebi.ac.uk/ena/browser/view/PRJEB43947

## DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

snapshot_date=${1:-'2021-11-29'}
dataset_name=${2:-'datahub_metadata'}
project_id=${3:-'prj-int-dev-covid19-nf-gls'}

output_dir="${DIR}/${snapshot_date}/output"
cd "${output_dir}" || exit

gsutil -m cp "${output_dir}/new_nanopore_metadata.tsv" "gs://${dataset_name}/nanopore_metadata.tsv"
gsutil -m cp "${output_dir}/all_illumina_metadata.tsv" "gs://${dataset_name}/illumina_metadata.tsv"

bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.nanopore_metadata" "gs://${dataset_name}/nanopore_metadata.tsv"

# There are bad collcetioon_date values. Use STRING instead of DATE.
bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.illumina_metadata" "gs://${dataset_name}/illumina_metadata.tsv"

#awk 'BEGIN{ FS=OFS="\t" }{print $0, "2021-11-29"}' /Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/2021-11-29/output/new_nanopore_metadata.tsv > /Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/2021-11-29/output/new_nanopore_metadata.1.tsv
#awk 'BEGIN{ FS=OFS="\t" }{print $0, "2021-11-29"}' /Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/2021-11-29/output/all_illumina_metadata.tsv > /Users/davidyuan/IdeaProjects/davidyuyuan/ENA_SARS_Cov2_nanopore/results/2021-11-29/output/all_illumina_metadata.1.tsv