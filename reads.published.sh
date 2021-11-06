#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22%20OR%20instrument_platform%3D%22ILLUMINA%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=100' "https://www.ebi.ac.uk/ena/portal/api/search" >all.tsv

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
  "https://www.ebi.ac.uk/ena/portal/api/search" >nanopore.index.tsv

curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d 'result=read_run&query=tax_tree(2697049)%20AND%20(instrument_platform%3D%22ILLUMINA%22)&fields=instrument_platform%2Cinstrument_model%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Ccountry%2Ccollection_date%2Cisolate%2Cstrain&format=tsv&limit=0' \
  "https://www.ebi.ac.uk/ena/portal/api/search" >illumina.index.tsv

export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA0MTE1fS5hZjhlZjYyMmI2NWFlZjE0YzY1ZWY5ODNiMjk4YmEwOWE4NzgzNzFi
export GOOGLE_APPLICATION_CREDENTIALS=~/gcp-nf/nextflow-service-account-private-key.json

~/gcp-nf/nextflow -C ~/gcp-nf/gls/nextflow.config run ~/ENA_SARS_Cov2_nanopore/nanopore_input.nf -w gs://prj-int-dev-covid19-nf-gls/nanopore_input -profile gls -with-tower --resume
~/gcp-nf/nextflow -C ~/gcp-nf/gls/nextflow.config run ~/ENA_SARS_Cov2_nanopore/nanopore_nextflow.nf -w gs://prj-int-dev-covid19-nf-gls/nanopore -profile gls -with-tower --resume --RUN_ID 1
