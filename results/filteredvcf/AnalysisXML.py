import hashlib
from lxml import etree


class AnalysisXML:
    def __init__(self, alias, run_accession, analysis_date,
                 analysis_file, analysis_title, analysis_description,
                 analysis_attributes, sample_accession):
        self.alias = alias
        self.project_accession = 'PRJEB45554'
        self.run_accession = run_accession
        self.analysis_date = analysis_date
        self.analysis_file = analysis_file
        self.analysis_title = analysis_title
        self.analysis_description = analysis_description
        self.analysis_attributes = analysis_attributes
        self.sample_accession = sample_accession
        self.centre_name = 'VEO'

    def build_xml(self):
        analysis_set = etree.Element('ANALYSIS_SET')
        analysis_xml = etree.ElementTree(analysis_set)

        analysis_elt = etree.SubElement(
            analysis_set, 'ANALYSIS', alias=self.alias,
            center_name=self.centre_name,
            analysis_date=self.analysis_date
        )
        # Title
        title_elt = etree.SubElement(analysis_elt, 'TITLE')
        title_elt.text = self.analysis_title

        # Description
        description_elt = etree.SubElement(analysis_elt, 'DESCRIPTION')
        description_elt.text = self.analysis_description

        # References on study, sample and run
        etree.SubElement(analysis_elt, 'STUDY_REF',
                         accession=self.project_accession)
        etree.SubElement(analysis_elt, 'SAMPLE_REF',
                         accession=self.sample_accession)
        etree.SubElement(analysis_elt, 'RUN_REF', accession=self.run_accession)

        # analysis type
        analyses_type_elt = etree.SubElement(analysis_elt, 'ANALYSIS_TYPE')
        etree.SubElement(analyses_type_elt, 'COVID19_FILTERED_VCF')

        # Files
        files_elt = etree.SubElement(analysis_elt, 'FILES')
        etree.SubElement(files_elt, 'FILE', filename=self.analysis_file,
                         filetype="vcf", checksum_method="MD5",
                         checksum=f"{self.get_checksum()}")

        # Analysis attributes
        analysis_attributes_elt = etree.SubElement(analysis_elt,
                                                   'ANALYSIS_ATTRIBUTES')
        for tag, value in self.analysis_attributes.items():
            analysis_attribute_elt = etree.SubElement(analysis_attributes_elt,
                                                      'ANALYSIS_ATTRIBUTE')
            etree.SubElement(analysis_attribute_elt, 'TAG').text = tag
            etree.SubElement(analysis_attribute_elt, 'VALUE').text = value

        analysis_xml.write(f"analysis_{self.analysis_date}.xml",
                           pretty_print=True,
                           xml_declaration=True, encoding='UTF-8')

    def get_checksum(self):
        return hashlib.md5(open(self.analysis_file, 'rb').read()).hexdigest()
