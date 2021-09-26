"""
functions to work with orcid data
Russ Poldrack, May 2020
"""

import pandas as pd


def get_dois_from_orcid_record(orcid_data):
    dois = []
    for g in orcid_data['activities-summary']['works']['group']:
        for p in g['work-summary']:
            doi = None
            for eid in p['external-ids']['external-id']:
                if eid['external-id-type'] == 'doi':
                    doi = eid['external-id-value'].replace('http://dx.doi.org/', '')
            if doi is not None:
                dois.append(doi.lower())
    return(list(set(dois)))


def get_orcid_education(orcid_data):
    education_df = pd.DataFrame(columns=['institution', 'degree', 'dept', 'city',
                                         'start_date', 'end_date'])
    for ctr, e in enumerate(orcid_data['activities-summary']['educations']['affiliation-group']):
        affiliation_info = parse_orcid_affiliation_record(
            e['summaries'][0]['education-summary'])
        education_df.loc[ctr, :] = affiliation_info

    for e in orcid_data['activities-summary']['qualifications']['affiliation-group']:
        ctr += 1
        affiliation_info = parse_orcid_affiliation_record(
            e['summaries'][0]['qualification-summary'])
        education_df.loc[ctr, :] = affiliation_info

    education_df = education_df.sort_values('start_date')
    return(education_df)


def parse_orcid_affiliation_record(record):
    institution = record['organization']['name']
    city = record['organization']['address']['city'] + ', ' + record['organization']['address']['region']
    start_date = record['start-date']['year']['value']
    end_date = record['end-date']['year']['value']
    degree = record['role-title']
    dept = record['department-name']
    return([institution, degree, dept, city, start_date, end_date])


def get_orcid_funding(orcid_data):
    funding_df = pd.DataFrame(columns=['organization', 'id', 'title', 'role',
                                       'start_date', 'end_date', 'url'])
    for ctr, e in enumerate(orcid_data['activities-summary']['fundings']['group']):
        funding_df.loc[ctr, :] = parse_orcid_funding_record(
            e['funding-summary'][0])

    return(funding_df)


def parse_orcid_funding_record(record):
    id = record['external-ids']['external-id'][0]['external-id-value']
    start_date = record['start-date']['year']['value']
    if record['end-date'] is not None:
        end_date = record['end-date']['year']['value']
    else:
        end_date = 'present'
    url = record['external-ids']['external-id'][0]['external-id-url']['value']
    return([record['organization']['name'], id,
           record['title']['title']['value'],
           '', start_date, end_date, url])


def get_orcid_employment(orcid_data):
    employment_df = pd.DataFrame(columns=['institution', 'role', 'dept', 'city',
                                          'start_date', 'end_date'])
    for ctr, e in enumerate(orcid_data['activities-summary']['employments']['affiliation-group']):
        employment_df.loc[ctr, :] = parse_orcid_employment_record(
            e['summaries'][0]['employment-summary'])
    employment_df = employment_df.sort_values('start_date', ascending=False)
    return(employment_df)


def parse_orcid_employment_record(record):
    institution = record['organization']['name']
    city = record['organization']['address']['city'] + ', ' + record['organization']['address']['region']
    start_date = record['start-date']['year']['value']
    if record['end-date'] is not None:
        end_date = record['end-date']['year']['value']
    else:
        end_date = 'present'
    role = record['role-title']
    dept = record['department-name']
    return([institution, role, dept, city, start_date, end_date])


def get_orcid_distinctions(orcid_data):
    distinctions_df = pd.DataFrame(columns=[
        'organization', 'title',
        'city', 'start_date', 'end_date',
        'distinction_type'])

    for ctr, e in enumerate(orcid_data['activities-summary']['distinctions']['affiliation-group']):
        distinctions_df.loc[ctr, :] = parse_orcid_distinctions_record(
            e['summaries'][0]['distinction-summary'],
            distinction_type='Honor')

    for e in orcid_data['activities-summary']['invited-positions']['affiliation-group']:
        ctr += 1
        distinctions_df.loc[ctr, :] = parse_orcid_distinctions_record(
            e['summaries'][0]['invited-position-summary'],
            distinction_type='Visiting position')

    distinctions_df = distinctions_df.sort_values('start_date', ascending=False)
    return(distinctions_df)


def parse_orcid_distinctions_record(record, distinction_type):
    organization = record['organization']['name']
    start_date = record['start-date']['year']['value']
    end_date = record['end-date']['year']['value'] if record['end-date'] is not None else ''
    role = record['role-title']
    city = record['organization']['address']['city']
    if record['organization']['address']['region'] is not None:
        city = city + ', ' + record['organization']['address']['region']
    return([organization, role, city, start_date, end_date, distinction_type])


def get_orcid_memberships(orcid_data):
    memberships_df = pd.DataFrame(columns=['organization'])

    for ctr, e in enumerate(orcid_data['activities-summary']['memberships']['affiliation-group']):
        s = e['summaries'][0]['membership-summary']
        organization = s['organization']['name']
        memberships_df.loc[ctr, :] = [organization]
    memberships_df = memberships_df.sort_values('organization')
    return(memberships_df)


def get_orcid_service(orcid_data):
    service_df = pd.DataFrame(columns=[
        'organization', 'role', 'start_date', 'end_date'])

    for ctr, e in enumerate(orcid_data['activities-summary']['services']['affiliation-group']):
        service_df.loc[ctr, :] = parse_orcid_service_record(e['summaries'][0]['service-summary'])

    service_df = service_df.sort_values('start_date', ascending=False)
    return(service_df)


def parse_orcid_service_record(record):
    organization = record['organization']['name']
    role = record['role-title']
    start_date = record['start-date']['year']['value']
    if record['end-date'] is not None:
        end_date = record['end-date']['year']['value']
    else:
        end_date = 'present'
    return([organization, role, start_date, end_date])
