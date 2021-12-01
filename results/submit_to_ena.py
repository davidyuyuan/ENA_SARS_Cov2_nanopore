import sys
from AnalysisXML import AnalysisXML
from SubmissionXML import SubmissionXML
from datetime import datetime
from lxml import etree
import subprocess

to_be_submitted = sys.argv[1]
snapshot_date = sys.argv[2]

timestamp = datetime.now()
analysis_date = timestamp.strftime("%Y-%m-%dT%H:%M:%S")
errors = open('errors.txt', 'w')
with open(to_be_submitted, 'r') as f:
    for line in f:
        line = line.rstrip()
        data = line.split("\t")
        analysis_file, run_accession, sample_accession, platform_name = data[1:]
        if platform_name == 'ILLUMINA':
            pipeline_name = 'COVID-19 Sequence Analysis Workflow'
        else:
            pipeline_name = 'Nanopore Analysis Workflow'
        alias = f"covid_sequence_analysis_workflow_{run_accession}_{analysis_date}"
        analysis_title = f"VEO SARS-CoV-2 systematically called variant data " \
                         f"of public run {run_accession} and sample " \
                         f"{sample_accession}, part of the snapshot " \
                         f"generated on 23/08/2021."
        analysis_description = f"All public SARS-CoV-2 INSDC raw read data " \
                               f"is streamed through the COVID Sequence " \
                               f"Analysis Workflow to produce a set of " \
                               f"uniform variant calls. For more information " \
                               f"on the pipeline, visit " \
                               f"https://github.com/enasequence/" \
                               f"covid-sequence-analysis-workflow. " \
                               f"VEO: https://www.veo-europe.eu/"
        analysis_attributes = {
            'PIPELINE_NAME': pipeline_name,
            'PIPELINE_VERSION': 'v1',
            'SUBMISSION_TOOL': 'ENA_Analysis_Submitter',
            'SUBMISSION_TOOL_VERSION': '1.0.0',
            'SNAPSHOT_DATE': '29/11/2021'}

        analysis_xml_obj = AnalysisXML(
            alias=alias, run_accession=run_accession,
            analysis_date=analysis_date, analysis_file=analysis_file,
            analysis_title=analysis_title,
            analysis_description=analysis_description,
            analysis_attributes=analysis_attributes,
            sample_accession=sample_accession)
        analysis_xml_obj.build_xml()

        submission_xml_obj = SubmissionXML(
            alias=alias, analysis_date=analysis_date)
        submission_xml_obj.build_xml()

        upload_process = subprocess.run(
            f"curl -T {analysis_file}  "
            f"ftp://webin.ebi.ac.uk --user Webin-56719:!-Content-!",
            shell=True, capture_output=True)
        checksum_process = subprocess.run(
            f"curl -s ftp://webin.ebi.ac.uk/{analysis_file} "
            f"--user Webin-56719:!-Content-! | md5sum | cut -f1 -d ' '",
            shell=True, capture_output=True)

        original_md5 = analysis_xml_obj.get_checksum()
        uploaded_md5 = checksum_process.stdout.decode('UTF-8').rstrip()

        if original_md5 == uploaded_md5:
            submit_process = subprocess.run(
                f'curl -u Webin-56719:!-Content-! '
                f'-F "SUBMISSION=@submission_{analysis_date}.xml" '
                f'-F "ANALYSIS=@analysis_{analysis_date}.xml" '
                f'"https://www.ebi.ac.uk/ena/submit/drop-box/submit/"',
                shell=True, capture_output=True)
            try:
                root = etree.fromstring(submit_process.stdout)
                if root.get('success') == 'true':
                    print(f"Successfully submitted: {run_accession}")
                else:
                    print(submit_process.stdout)
                    print(f"submission error: {run_accession}")
                    errors.write(f"{run_accession}\n")
            except etree.XMLSyntaxError:
                print(f"submission error: {run_accession}")
                errors.write(f"{run_accession}\n")
        else:
            print(f"md5 error: {run_accession}")
            errors.write(f"{run_accession}\n")
errors.close()
