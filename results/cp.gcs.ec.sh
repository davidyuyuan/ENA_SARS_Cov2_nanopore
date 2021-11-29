#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Install Minoconda with Python 3.9
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

git init

# Install Google Cloud SSDK
sudo tee -a /etc/yum.repos.d/google-cloud-sdk.repo << EOM
[google-cloud-sdk]
name=Google Cloud SDK
baseurl=https://packages.cloud.google.com/yum/repos/cloud-sdk-el8-x86_64
enabled=1
gpgcheck=1
repo_gpgcheck=0
gpgkey=https://packages.cloud.google.com/yum/doc/yum-key.gpg
       https://packages.cloud.google.com/yum/doc/rpm-package-key.gpg
EOM
sudo dnf install google-cloud-sdk
gcloud init

# Attach and mount Cinder dy-nextflow: 8 TiB SSD
lsblk
sudo mkfs.ext4 -m 0 -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/vdb
sudo mkdir -p /mnt/result
sudo mount -o discard,defaults /dev/vdb /mnt/result
sudo chmod a+w /mnt/result
sudo blkid /dev/vdb
sudo vi /etc/fstab
#UUID=5564fd53-d902-4124-8187-70c25c469251 /mnt/result ext4 discard,defaults 0 2

# Pull new results and duplicated results separately
nohhup gsutil -m cp -n -I /mnt/result/from_gcs/nanopore/new < /mnt/result/from_gcs/nanopore/new_results.txt
nohhup gsutil -m cp -n -I /mnt/result/from_gcs/nanopore/duplicated < /mnt/result/from_gcs/nanopore/duplicated_results.txt


##!/usr/bin/env bash
#
#export AWS_ACCESS_KEY_ID=D0FVX5OMT6UWFNN5DVFA
#export AWS_SECRET_ACCESS_KEY=GwIROUcKvad2bjhuEbVu6TNQr7kPToeisrPABk3Z
#export AWS_DEFAULT_REGION=us-east-1
#
## DIR where the current script resides
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#
#while IFS="" read -r p || [ -n "$p" ]
#do
##   printf '%s\n' "$p"
##   printf '%s\n' $(basename "$p")
#    gsutil -m cp "$p" "${DIR}"
#    aws s3 --endpoint-url https://uk1s3.embassy.ebi.ac.uk cp $(basename "$p") s3://results_from_gcs
#    rm $(basename "$p")
#    # ls -ltr
#    # aws s3 --endpoint-url https://uk1s3.embassy.ebi.ac.uk ls s3://results_from_gcs
#done < "${DIR}/short_results.txt"
