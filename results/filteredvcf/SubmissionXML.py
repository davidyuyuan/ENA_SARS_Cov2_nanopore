from lxml import etree


class SubmissionXML:
    def __init__(self, alias, analysis_date):
        self.alias = alias
        self.action = 'ADD'
        self.centre_name = 'VEO'
        self.analysis_date = analysis_date
        self.source_xml = f'analysis_{self.analysis_date}.xml'
        self.schema = 'analysis'

    def build_xml(self):
        submission_set_elt = etree.Element('SUBMISSION_SET')
        submission_xml = etree.ElementTree(submission_set_elt)

        submission_elt = etree.SubElement(submission_set_elt, 'SUBMISSION',
                                          alias=self.alias,
                                          center_name=self.centre_name)
        actions_elt = etree.SubElement(submission_elt, 'ACTIONS')
        action_elt = etree.SubElement(actions_elt, 'ACTION')
        etree.SubElement(action_elt, self.action, source=self.source_xml,
                         schema=self.schema)
        submission_xml.write(f"submission_{self.analysis_date}.xml",
                             pretty_print=True, xml_declaration=True,
                             encoding='UTF-8')
