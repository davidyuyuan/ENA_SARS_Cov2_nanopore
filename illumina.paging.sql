SELECT * FROM `prj-int-dev-covid19-nf-gls.datahub_metadata.illumina_index` T1
WHERE T1.fastq_ftp IS NOT NULL AND T1.run_accession NOT IN
    (SELECT T2.run_id FROM `prj-int-dev-covid19-nf-gls.datahub_metadata.illumina_metadata` T2)

SELECT * FROM `prj-int-dev-covid19-nf-gls.datahub_metadata.nanopore_index` T1
WHERE T1.fastq_ftp IS NOT NULL AND T1.run_accession NOT IN
    (SELECT T2.run_id FROM `prj-int-dev-covid19-nf-gls.datahub_metadata.nanopore_metadata` T2)
