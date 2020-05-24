#!/usr/bin/env python
"""
generate a LaTeX CV using PubMed and ORCID along with other resources

Russ Poldrack, May 2020
"""

import os
import shutil
import pkg_resources
from autocv.autocv import get_params, get_orcid_data, get_orcid_education,\
    get_orcid_employment, get_orcid_distinctions, get_orcid_memberships,\
    get_orcid_service, get_orcid_dois
from autocv.autocv import get_funding_from_csv, get_teaching, get_patents,\
    drop_excluded_pubs, add_additional_pubs_from_csv, get_links
from autocv.autocv import render_education, render_distinctions,\
    render_employment, render_memberships,\
    render_service, render_funding, render_teaching, \
    render_editorial, render_heading, render_patents, render_pubs
from autocv.autocv import get_google_scholar_record
from autocv.autocv import get_pubmed_records, get_pubmed_pubs
from autocv.autocv import get_crossref_records, process_crossref_records
from autocv.autocv import write_pubs, write_presentations, write_talks


if __name__ == "__main__":

    params = get_params()

    # get orcid data
    print('reading data for ORCID id:', params['orcid'])
    orcid_data = get_orcid_data(params['orcid'])

    education_df = get_orcid_education(orcid_data)
    render_education(education_df)

    employment_df = get_orcid_employment(orcid_data)
    render_employment(employment_df)

    distinctions_df = get_orcid_distinctions(orcid_data)
    render_distinctions(distinctions_df)

    memberships_df = get_orcid_memberships(orcid_data)
    render_memberships(memberships_df)

    service_df = get_orcid_service(orcid_data)
    render_service(service_df)

    funding_df = get_funding_from_csv()
    render_funding(funding_df)

    teaching_dict = get_teaching()
    render_teaching(teaching_dict)

    render_editorial()

    render_heading(params)

    patent_list = get_patents(params['lastname'], params['firstname'])
    if len(patent_list) > 0:
        render_patents(patent_list)

    gscholar = get_google_scholar_record(params['firstname'], params['lastname'])
    download_pubs = True  # use for testing

    pubmed_records = get_pubmed_records(params)
    pubmed_pubs = get_pubmed_pubs(pubmed_records)
    orcid_dois = get_orcid_dois(orcid_data)

    all_dois = list(set(orcid_dois + list(pubmed_pubs.keys())))
    print('found %d total records across Pubmed and ORCID' % len(all_dois))

    if download_pubs:
        crossref_records = get_crossref_records(all_dois)

    # use pubmed info if it exists
    pubs = process_crossref_records(crossref_records, pubmed_pubs)

    # drop excluded pubs
    pubs = drop_excluded_pubs(pubs)

    # load missing pubs from file
    pubs = add_additional_pubs_from_csv(pubs)

    # get additional links and add to pubs
    links = get_links('links.csv')

    for linktype in links:
        for id in links[linktype]:
            if 'links' not in pubs[id]:
                pubs[id]['links'] = {}
            pubs[id]['links'][linktype] = links[linktype][id]

    latex_lines = render_pubs(pubs, gscholar)
    write_pubs(latex_lines)
    write_presentations()
    write_talks()

    # if template is not in current directory then copy it here
    if not os.path.exists('autocv_template.tex'):
        tempfile = pkg_resources.resource_filename('autocv', 'templates/autocv_template.tex')
        shutil.copy(tempfile, '.')
