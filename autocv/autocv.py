"""
generate a LaTeX CV using PubMed and ORCID along with other resources

Russ Poldrack, May 2020
"""

import pandas as pd
import numpy as np
from Bio import Entrez
import json
import os
from pprint import pprint
import random
import string
import requests
from crossref.restful import Works
import pypatent
import scholarly
from collections import OrderedDict


# load parameters
def get_params(params_file='params.json'):
    if os.path.exists(params_file):
        with open(params_file) as f:
            params = json.load(f)
    else:
        raise FileNotFoundError('Please create a json file called params.json containing the fields email (with your email address), orcid (with your ORCID id) and query (with your pubmed query)- see documentation for help')
    required_fields = ['address', 'lastname', 'firstname', 'email', 'orcid', 'query', 'url', 'phone']
    for field in required_fields:
        assert field in params
    return(params)


def get_orcid_data(id):
    resp = requests.get("http://pub.orcid.org/%s" % id,
                        headers={'Accept': 'application/orcid+json'})
    orcid_data = resp.json()
    return(orcid_data)


def get_orcid_education(orcid_data):
    education_df = pd.DataFrame(columns=['institution', 'degree', 'dept', 'city',
                                         'start_date', 'end_date'])
    ctr = 0
    for e in orcid_data['activities-summary']['educations']['affiliation-group']:
        s = e['summaries'][0]['education-summary']
        institution = s['organization']['name']
        city = s['organization']['address']['city'] + ', ' + s['organization']['address']['region']
        start_date = s['start-date']['year']['value']
        end_date = s['end-date']['year']['value']
        degree = s['role-title']
        dept = s['department-name']
        education_df.loc[ctr, :] = [institution, degree, dept, city, start_date, end_date]
        ctr += 1

    for e in orcid_data['activities-summary']['qualifications']['affiliation-group']:
        s = e['summaries'][0]['qualification-summary']
        institution = s['organization']['name']
        city = s['organization']['address']['city'] + ', ' + s['organization']['address']['region']
        start_date = s['start-date']['year']['value']
        end_date = s['end-date']['year']['value']
        degree = s['role-title']
        dept = s['department-name']
        education_df.loc[ctr, :] = [institution, degree, dept, city, start_date, end_date]
        ctr += 1

    education_df = education_df.sort_values('start_date')
    return(education_df)


def get_orcid_funding(orcid_data):
    funding_df = pd.DataFrame(columns=['organization', 'id', 'title', 'role',
                                       'start_date', 'end_date', 'url'])
    ctr = 0
    for e in orcid_data['activities-summary']['fundings']['group']:
        s = e['funding-summary'][0]
        id = s['external-ids']['external-id'][0]['external-id-value']
        start_date = s['start-date']['year']['value']
        if s['end-date'] is not None:
            end_date = s['end-date']['year']['value']
        else:
            end_date = 'present'
        url = s['external-ids']['external-id'][0]['external-id-url']['value']
        funding_df.loc[ctr, :] = [s['organization']['name'], id,
                                  s['title']['title']['value'],
                                  '',
                                  start_date,
                                  end_date, url]
        ctr += 1
        return(funding_df)


def get_orcid_employment(orcid_data):
    employment_df = pd.DataFrame(columns=['institution', 'role', 'dept', 'city',
                                          'start_date', 'end_date'])
    ctr = 0
    for e in orcid_data['activities-summary']['employments']['affiliation-group']:
        s = e['summaries'][0]['employment-summary']
        institution = s['organization']['name']
        city = s['organization']['address']['city'] + ', ' + s['organization']['address']['region']
        start_date = s['start-date']['year']['value']
        if s['end-date'] is not None:
            end_date = s['end-date']['year']['value']
        else:
            end_date = 'present'
        role = s['role-title']
        dept = s['department-name']
        employment_df.loc[ctr, :] = [institution, role, dept, city, start_date, end_date]
        ctr += 1
    employment_df = employment_df.sort_values('start_date', ascending=False)
    return(employment_df)


def get_orcid_distinctions(orcid_data):
    distinctions_df = pd.DataFrame(columns=['organization', 'title', 'city', 'start_date', 'end_date'])

    for ctr, e in enumerate(orcid_data['activities-summary']['distinctions']['affiliation-group']):
        s = e['summaries'][0]['distinction-summary']
        organization = s['organization']['name']
        start_date = s['start-date']['year']['value']
        if s['end-date'] is not None:
            end_date = s['end-date']['year']['value']
        else:
            end_date = ''
        role = s['role-title']
        distinctions_df.loc[ctr, :] = [organization, role, '', start_date, end_date]

    for e in orcid_data['activities-summary']['invited-positions']['affiliation-group']:
        s = e['summaries'][0]['invited-position-summary']
        institution = s['organization']['name']
        city = s['organization']['address']['city']
        if s['organization']['address']['region'] is not None:
            city = city + ', ' + s['organization']['address']['region']
        start_date = s['start-date']['year']['value']
        if s['end-date'] is not None:
            end_date = s['end-date']['year']['value']
        else:
            end_date = 'present'
        role = s['role-title']
        distinctions_df.loc[ctr, :] = [institution, role, city, start_date, end_date]
        ctr += 1

    distinctions_df = distinctions_df.sort_values('start_date', ascending=False)
    return(distinctions_df)


def get_orcid_memberships(orcid_data):
    memberships_df = pd.DataFrame(columns=['organization'])

    for ctr, e in enumerate(orcid_data['activities-summary']['memberships']['affiliation-group']):
        s = e['summaries'][0]['membership-summary']
        organization = s['organization']['name']
        memberships_df.loc[ctr, :] = [organization]
    memberships_df = memberships_df.sort_values('organization')
    return(memberships_df)


def get_orcid_service(orcid_data):
    service_df = pd.DataFrame(columns=['organization'])

    for ctr, e in enumerate(orcid_data['activities-summary']['services']['affiliation-group']):
        s = e['summaries'][0]['service-summary']
        service_df.loc[ctr, 'organization'] = s['organization']['name']
        service_df.loc[ctr, 'start_date'] = s['start-date']['year']['value']
        if s['end-date'] is not None:
            service_df.loc[ctr, 'end_date'] = s['end-date']['year']['value']
        else:
            service_df.loc[ctr, 'end_date'] = 'present'
        service_df.loc[ctr, 'role'] = s['role-title']

    service_df = service_df.sort_values('start_date', ascending=False)
    return(service_df)


# get DOIs for pubs in orcid
def get_orcid_dois(orcid_data):
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


# get records for pubs in PubMed
def get_pubmed_records(params):
    Entrez.email = params['email']
    print(f'using {params["email"]} for Entrez service')
    query = params['query']
    print('searching for', query)
    retmax = 1000
    handle = Entrez.esearch(db="pubmed", retmax=retmax, term=query)
    record = Entrez.read(handle)
    handle.close()
    pmids = [int(i) for i in record['IdList']]
    print('found %d matches' % len(pmids))

    # load full records
    handle = Entrez.efetch(db="pubmed", id=",".join(['%d' % i for i in pmids]),
                           retmax=retmax, retmode="xml")
    records = Entrez.read(handle)
    print('found %d full pubmed records' % len(records['PubmedArticle']))
    return(records)


def get_pubmed_pubs(pubmed_records):
    pubs = {}
    for i in pubmed_records['PubmedArticle']:
        pmc = None
        doi = None
        pmid = int(i['MedlineCitation']['PMID'])
        for j in i['PubmedData']['ArticleIdList']:
            if j.attributes['IdType'] == 'doi':
                doi = str(j).lower().replace('http://dx.doi.org/', '')
            if j.attributes['IdType'] == 'pmc':
                pmc = str(j)

        if doi is None:
            continue

        pubs[doi] = {'pmid': pmid, 'PMC': pmc, 'type': 'journal-article'}
        # get some other useful stuff while we are here
        # pubmed seems to have better info than crossref
        pubs[doi]['journal'] = i['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        if 'Year' in i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            pubs[doi]['year'] = int(i['MedlineCitation']['Article'][
                'Journal']['JournalIssue']['PubDate']['Year'])
        elif 'MedlineDate' in i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            pubs[doi]['year'] = int(i['MedlineCitation']['Article'][
                'Journal']['JournalIssue']['PubDate']['MedlineDate'].split(' ')[0])
        if 'Volume' in i['MedlineCitation']['Article']['Journal']['JournalIssue']:
            pubs[doi]['volume'] = i['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
        pubs[doi]['title'] = i['MedlineCitation']['Article']['ArticleTitle']

        if 'Pagination' in i['MedlineCitation']['Article']:
            pubs[doi]['page'] = i['MedlineCitation']['Article']['Pagination']['MedlinePgn']
        authorlist = []
        if 'AuthorList' in i['MedlineCitation']['Article']:
            for author in i['MedlineCitation']['Article']['AuthorList']:
                if 'LastName' in author and 'Initials' in author:
                    authorlist.append(' '.join([author['LastName'], author['Initials']]))
            pubs[doi]['authors'] = ', '.join(authorlist)

    return(pubs)


def get_crossref_records(dois):
    works = Works()
    crossref_records = {}
    print('searching crossref for all DOIs, this might take a few minutes...')
    for doi in dois:
        r = works.doi(doi)
        if r is not None:
            crossref_records[doi] = r
        else:
            print('missing crossref record for', doi)
    return(crossref_records)


def process_crossref_records(crossref_records, pubs,
                             etal_thresh=10, exclude_preprints=True,
                             exclude_books=True, exclude_translations=True):
    # use DOI as dictionary keys
    for r in crossref_records:
        if exclude_preprints and crossref_records[r]['type'] == 'posted-content':
            continue
        # books seem to be goofy with crossref
        if exclude_books and crossref_records[r]['type'] == 'book':
            continue
        if 'author' not in crossref_records[r]:  # can happen for errata
            print('no author for ', r)
            continue
        if 'translator' in crossref_records[r]:
            continue

        if isinstance(crossref_records[r]['title'], list):
            crossref_records[r]['title'] = crossref_records[r]['title'][0]
        if crossref_records[r]['title'].find('Corrigend') > -1:
            continue

        if r not in pubs:
            pubs[r] = {}

        # don't replace pubmed info if it already exists
        for field in ['title', 'volume', 'page', 'type', 'publisher']:
            if field in crossref_records[r] and field not in pubs[r]:
                f = crossref_records[r][field]
                if isinstance(f, list):
                    f = f[0]
                pubs[r][field] = f

        # get the title
        if 'journal' not in pubs[r] and len(crossref_records[r]['container-title']) > 0:
            pubs[r]['journal'] = crossref_records[r]['container-title'][0]

        # date can show up in two different places!
        if 'year' not in pubs[r]:
            if 'published-print' in crossref_records[r]:
                year = crossref_records[r]['published-print']['date-parts'][0][0]
            elif 'journal-issue' in crossref_records[r]:
                journal_issue = crossref_records[r]['journal-issue']
                if 'published-print' in journal_issue:
                    year = journal_issue['published-print']['date-parts'][0][0]
                else:
                    year = journal_issue['published-online']['date-parts'][0][0]

            pubs[r]['year'] = int(year)

        # convert author list to pubmed format
        if 'authors' in pubs[r]:
            authors = [i.lstrip(' ') for i in pubs[r]['authors'].split(',')]
        else:
            authors = []
            for author in crossref_records[r]['author']:
                if 'given' not in author or 'family' not in author:
                    continue
                given_split = author['given'].split(' ')
                if len(given_split) > 1:
                    initials = ''.join([i[0] for i in given_split])
                else:
                    initials = given_split[0][0]
                entry = '%s %s' % (author['family'], initials)
                authors.append(entry)
        if len(authors) > etal_thresh:
            pubs[r]['authors'] = '%s et al.' % authors[0]
        else:
            pubs[r]['authors'] = ', '.join(authors)
        pubs[r]['idType'] = 'DOI'
    return(pubs)


def get_google_scholar_record(firstname, lastname):
    search_query = scholarly.scholarly.search_author(' '.join([firstname, lastname]))
    author = next(search_query).fill()
    return(author)


# get code/data/osf links
def get_links(link_file, verbose=False):
    links = {}
    if os.path.exists(link_file):
        data_df = pd.read_csv(link_file)
        for i in data_df.index:
            linktype = data_df.loc[i, 'type']
            id = data_df.loc[i, 'DOI']
            if linktype not in links:
                links[linktype] = {}
            links[linktype][id] = data_df.loc[i, 'url']
        if verbose:
            print('Links:')
            pprint(links)
    return(links)


def add_additional_pubs_from_csv(pubs, pubfile='additional_pubs.csv'):
    if os.path.exists(pubfile):
        addpubs = pd.read_csv(pubfile)
        addpubs = addpubs.fillna('')
        # resolve duplicate ISBNs (e.g. multiple chapters in a book)
        if 'ISBN' in addpubs.columns:
            isbn_loc = np.where(addpubs.columns == 'ISBN')[0][0]
            addpubs.loc[:, 'ISBN'] = [i.strip(' ') for i in addpubs.ISBN]
            for i in range(1, addpubs.shape[0]):
                if addpubs.iloc[i, isbn_loc] == '':
                    continue
                if addpubs.iloc[i, isbn_loc] in addpubs.iloc[:(i - 1), isbn_loc].tolist():
                    print('found match')
                    addpubs.iloc[i, isbn_loc] = addpubs.iloc[i, isbn_loc] + '-' + ''.join(
                        random.choice(string.ascii_lowercase) for i in range(3))
        for i in addpubs.index:
            # make a random string to stand in for pmid
            if addpubs.loc[i, 'DOI'] != '':
                id = addpubs.loc[i, 'DOI']
                idType = 'DOI'
            elif addpubs.loc[i, 'ISBN'] != '':
                id = addpubs.loc[i, 'ISBN']
                idType = 'ISBN'
            else:
                id = ''.join(random.choice(
                    string.ascii_lowercase) for i in range(8))
                idType = 'randomID'
            if id in pubs and pubs[id]['title'] == addpubs.loc[i, 'title']:
                print('found duplicate pub - skipping:', id)
                continue
            elif id in pubs:
                print('found duplicate id - modifying:', id)
                print(pubs[id])
                print('')
                id += '-' + ''.join(random.choice(
                    string.ascii_lowercase) for i in range(8))
            else:
                print('found unique id:', id)
                print('')
            pubs[id] = {idType: id}
            for c in addpubs.columns:
                if c in ['DOI', 'ISBN']:
                    continue
                entry = addpubs.loc[i, c]
                if isinstance(entry, str):
                    entry = entry.strip(' ')
                pubs[id][c] = entry
    return(pubs)


def add_books_from_csv(pubs, bookfile='books.csv'):
    if os.path.exists(bookfile):
        books = pd.read_csv(bookfile)
        books = books.fillna('')
        for i in books.index:
            # make a random string to stand in for pmid
            id = books.loc[i, 'ISBN']
            pubs[id] = {'type': 'book'}
            for c in books.columns:
                pubs[id][c] = books.loc[i, c]
    return(pubs)


def get_teaching(teaching_file='teaching.csv'):
    teaching_dict = OrderedDict()
    if os.path.exists(teaching_file):
        teaching_df = pd.read_csv(teaching_file)
        for i in teaching_df.index:
            coursetype = teaching_df.loc[i, 'type']
            if coursetype not in teaching_dict:
                teaching_dict[coursetype] = []
            teaching_dict[coursetype].append(teaching_df.loc[i, 'name'])

    return(teaching_dict)


def get_patents(lastname, firstname):
    results = pypatent.Search(lastname).as_list()
    mypatents = []
    for r in results:
        for i in r['inventors']:
            fn = i[0].split(' ')[0].lower()
            ln = i[1].lower()
            if fn == firstname.lower() and ln == lastname.lower():
                mypatents.append(r)
    return(mypatents)


def get_funding_from_csv(funding_file='funding.csv'):
    if os.path.exists(funding_file):
        grants_df = pd.read_csv(funding_file)
        return(grants_df)
    else:
        return(None)


def drop_excluded_pubs(pubs, exclusions_file='exclusions.txt'):
    if os.path.exists(exclusions_file):
        e = pd.read_csv(exclusions_file)
        for i in e.index:
            doi = e.loc[i, 'DOI']
            if doi in pubs:
                print('dropping excluded doi:', doi)
                del pubs[doi]
    return(pubs)


def render_patents(patent_list):
    with open('patents.tex', 'w') as f:
        f.write('\\section*{Patents}\n\\noindent\n\n')
        for p in patent_list:
            authorlist = []
            for a in p['inventors']:
                initials = [i[0] for i in a[0].split(' ')]
                ln = a[1]
                authorlist.append('%s %s' % (ln, ''.join(initials)))
            authors = ', '.join(authorlist)

            f.write('%s (%s) \\textit{%s} US Patent \# \\href{%s}{%s} \\vspace{2mm}\n\n' % ( # noqa
                authors,
                p['patent_date'],
                p['title'],
                p['url'],
                p['patent_num']
            ))


def render_education(education_df):
    with open('education.tex', 'w') as f:
        f.write('\\section*{Education and training}\n\\noindent\n\n')
        for i in education_df.index:

            f.write('\\textit{%s-%s}: %s (%s), %s, %s\n\n' % (
                education_df.loc[i, 'start_date'],
                education_df.loc[i, 'end_date'],
                education_df.loc[i, 'degree'],
                education_df.loc[i, 'dept'],
                education_df.loc[i, 'institution'],
                education_df.loc[i, 'city'],
            ))


def render_employment(employment_df):
    with open('employment.tex', 'w') as f:
        f.write('\\section*{Employment and professional affiliations}\n\\noindent\n\n')
        for i in employment_df.index:
            if employment_df.loc[i, 'dept'] is None:
                dept = ''
            else:
                dept = ' (%s)' % employment_df.loc[i, 'dept']

            f.write('\\textit{%s-%s}: %s%s, %s\n\n' % (
                employment_df.loc[i, 'start_date'],
                employment_df.loc[i, 'end_date'],
                employment_df.loc[i, 'role'],
                dept,
                employment_df.loc[i, 'institution'],
            ))


def render_distinctions(distinctions_df):
    with open('distinctions.tex', 'w') as f:
        f.write('\\section*{Honors and Awards}\n\\noindent\n\n')
        for i in distinctions_df.index:
            f.write('\\textit{%s}: %s, %s\n\n' % (
                distinctions_df.loc[i, 'start_date'],
                distinctions_df.loc[i, 'title'],
                distinctions_df.loc[i, 'organization'],
            ))


def render_service(service_df):
    with open('service.tex', 'w') as f:
        f.write('\\section*{Service}\n\\noindent\n\n')
        for i in service_df.index:
            f.write('%s, %s, %s-%s \n\n' % (
                service_df.loc[i, 'role'],
                service_df.loc[i, 'organization'],
                service_df.loc[i, 'start_date'],
                service_df.loc[i, 'end_date'],
            ))


def render_memberships(memberships_df):
    with open('memberships.tex', 'w') as f:
        f.write('\\section*{Professional societies}\n\\noindent\n\n')
        memberships = ', '.join(memberships_df.organization)
        f.write('%s\n\n' % memberships)


def render_teaching(teaching_dict):
    with open('teaching.tex', 'w') as f:
        f.write('\\section*{Teaching}\n\\noindent\n\n')
        for i in teaching_dict:
            f.write('\\textit{%s}: %s \\vspace{2mm}\n\n' % (i, ', '.join(teaching_dict[i])))


def render_heading(params):
    lines = ['\\reversemarginpar \n']
    lines.append('{\LARGE %s %s. %s}\\\\[4mm] \n' % (  # noqa
        params['firstname'].title(),
        params['middlename'].title()[0],
        params['lastname'].title()))
    lines.append('\\vspace{-1cm} \n\n')

    lines.append('\\begin{multicols}{2} \n')
    for a in params['address']:
        lines.append(a + '\\\\\n')
    lines.append('\\columnbreak \n\n')
    lines.append('Phone: %s \\\\\n' % params['phone'])
    lines.append('email: %s \\\\\n' % params['email'])
    lines.append('url: \\href{%s}{%s} \\\\\n' % (
        params['url'].replace('http://', ''), params['url']))
    if 'github' in params:
        lines.append('url: \\href{%s}{%s} \\\\\n' % (
            params['github'], params['github'].replace('http://', '')))
    if 'twitter' in params:
        lines.append('Twitter: %s \\\\\n' % params['twitter'])
    lines.append('ORCID: \\href{https://orcid.org/%s}{%s} \\\\\n' % (
        params['orcid'], params['orcid']))
    lines.append('\end{multicols}\n\n')  # noqa
    with open('header.tex', 'w') as f:
        for l in lines:
            f.write(l)


def render_editorial(editorial_file='editorial.csv'):
    if os.path.exists(editorial_file):
        editorial_df = pd.read_csv(editorial_file)
    else:
        return

    editorial_df = editorial_df.fillna('')
    editorial_dict = OrderedDict()
    for i in editorial_df.index:
        role = editorial_df.loc[i, 'role']
        if role not in editorial_dict:
            editorial_dict[role] = []
        if editorial_df.loc[i, 'dates'] != '':
            date_string = ' (%s)' % editorial_df.loc[i, 'dates']
        else:
            date_string = ''
        editorial_dict[role].append(editorial_df.loc[i, 'journal'].strip(' ') + date_string)

    with open('editorial.tex', 'w') as f:
        f.write('\\section*{Editorial Duties and Reviewing} \n\\noindent \n\n')
        for i in editorial_dict:
            f.write('\\textit{%s}: %s \n\n' % (i.strip(' '),
                    ', '.join(editorial_dict[i])))


def make_funding_line(funding_df, i, abbreviate=True):
    if funding_df.loc[i, 'organization'].find('National') == 0 and abbreviate:
        # abbreviate
        org_split = [x[0] for x in funding_df.loc[
            i, 'organization'].split(' ') if x not in ['of', 'for', 'and', 'on']]
        org = ''.join(org_split)
    else:
        org = funding_df.loc[i, 'organization']
    if funding_df.loc[i, 'id'] == '':
        idtext = ''
    else:
        if funding_df.loc[i, 'url'] != '':
            idtext = ' (\\href{%s}{%s})' % (
                funding_df.loc[i, 'url'],
                funding_df.loc[i, 'id'])
        else:
            idtext = ' (%s)' % funding_df.loc[i, 'id']
    line = '%s, %s%s, \\textit{%s}, %s-%s' % (
        funding_df.loc[i, 'role'],
        org,
        idtext,
        funding_df.loc[i, 'title'].title().strip(' '),
        funding_df.loc[i, 'start_date'],
        funding_df.loc[i, 'end_date']
    )
    return(line)


def render_funding(funding_df, abbreviate=True):
    with open('funding.tex', 'w') as f:
        f.write('\\section*{Research funding}\n\\noindent\n\n')
        funding_df = funding_df.fillna('')
        # first print active grants
        f.write('\\subsection*{Active:}\n\n')

        # first print PI
        active_pi_funding = funding_df.query('active == True and role == "Principal Investigator"')
        for i in active_pi_funding.index:
            line = make_funding_line(active_pi_funding, i, abbreviate)
            f.write('%s \\vspace{2mm}\n\n' % line)

        # first print PI
        active_coi_funding = funding_df.query('active == True and role != "Principal Investigator"')
        for i in active_coi_funding.index:
            line = make_funding_line(active_coi_funding, i, abbreviate)
            f.write('%s \\vspace{2mm}\n\n' % line)

        # then print completed grants
        completed_funding = funding_df.query('active == False')
        f.write('\\subsection*{Completed:}\n\n')
        for i in completed_funding.index:
            line = make_funding_line(completed_funding, i, abbreviate)
            f.write('%s \\vspace{2mm}\n\n' % line)


def get_pubs_by_year(pubs, year):
    year_pubs = {}
    for p in pubs:
        if pubs[p]['year'] == year:
            year_pubs[p] = pubs[p]
    return(year_pubs)


def make_article_reference(pub):
    line = pub['authors'] +\
        ' (%s). ' % pub['year'] +\
        pub['title'] +\
        ' \\textit{%s' % pub['journal']
    if 'volume' in pub:
        line += ', %s}' % pub['volume']
    else:
        line += '}'
    if 'page' in pub:
        if len(pub['page']) > 0:
            line += ', %s' % pub['page']
    line += '.'
    return(line)


def make_chapter_reference(pub):
    page_string = ''
    if 'page' in pub:
        if len(pub['page']) > 0:
            page_string = '(p. %s). ' % pub['page']
    line = pub['authors'] +\
        ' (%s). ' % pub['year'] +\
        pub['title'] +\
        '. In \\textit{%s.} %s%s.' % (
            pub['journal'],
            page_string,
            pub['publisher'].strip(' '))
    return(line)


def make_book_reference(pub):
    line = pub['authors'] +\
        ' (%s). ' % pub['year'] +\
        ' \\textit{%s}. ' % pub['title'].strip(' ').strip('.') + \
        pub['publisher'].strip(' ')
    line += '.'
    return(line)


def render_pubs(pubs, gscholar):
    # convert to latex

    # get list of years
    years = list(set([pubs[i]['year'] for i in pubs]))
    years.sort(reverse=True)

    latex_lines = ['\\section*{Publications (Google Scholar H-index = %d)}' % gscholar.hindex]

    for year in years:
        latex_lines.append('\\subsection*{%s}' % year)
        year_pubs = get_pubs_by_year(pubs, year)
        # first get alphabetical order of PMIDs
        pmid_df = pd.DataFrame({'author': ''}, index=list(year_pubs.keys()))
        for pmid in year_pubs:
            if 'authors' in year_pubs[pmid]:
                pmid_df.loc[pmid, 'author'] = year_pubs[pmid]['authors']
            else:
                print('missing author:', pmid)
        pmid_df.sort_values('author', inplace=True)
        pmids_sorted = list(pmid_df.index)

        # get the reference line for each
        # NOTE: "pmid" here is a misnomer it's actually the DOI
        for pmid in pmids_sorted:
            for field in pubs[pmid]:
                if hasattr(pubs[pmid][field], 'replace'):
                    pubs[pmid][field] = pubs[pmid][field].replace(' &', ' \&')  # noqa
            pubs[pmid]['DOI'] = pmid

            line = None
            if pubs[pmid]['type'] in ['book', 'monograph']:
                line = make_book_reference(pubs[pmid])
            elif pubs[pmid]['type'] == 'book-chapter':
                line = make_chapter_reference(pubs[pmid])
            elif pubs[pmid]['type'] in ['proceedings-article', 'journal-article']:
                line = make_article_reference(pubs[pmid])
            if line is None:
                continue

            if 'PMC' in pubs[pmid]:
                if pubs[pmid]['PMC'] is not None:
                    line += ' \\href{https://www.ncbi.nlm.nih.gov/pmc/articles/%s}{OA}' % pubs[pmid]['PMC']

            if 'links' in pubs[pmid]:
                for linktype in pubs[pmid]['links']:
                    line += ' \\href{%s}{%s}' % (pubs[pmid]['links'][linktype],
                                                 linktype)

            # TBD: need to filter out non-DOIs
            if pubs[pmid]['type'] not in ['book', 'monograph'] and 'DOI' in pubs[pmid]:
                line += ' \\href{http://dx.doi.org/%s}{DOI}' % pubs[pmid]['DOI']

            line += ' \\vspace{2mm}'
            latex_lines.append(line)

    return(latex_lines)


def write_pubs(latex_lines):
    print('writing pubs to pubs.tex')
    with open('pubs.tex', 'w') as f:
        for l in latex_lines:
            f.write(l.replace('_', '\_') + '\n\n') # noqa


def write_presentations(presentations_file='conference.csv'):
    if os.path.exists(presentations_file):
        presentations = pd.read_csv(presentations_file, index_col=0)
        presentations = presentations.sort_values('year', ascending=False)

        print('writing presentations to presentations.tex')
        with open('presentations.tex', 'w') as f:
            f.write('\\section*{Conference Presentations}\n\\noindent\n\n')
            for i in presentations.index:
                entry = presentations.loc[i, :]
                title = entry.title.strip('.')
                location = entry.location.strip(' ').strip('.')
                term = '\\vspace{2mm}'
                line = '%s (%s). \\emph{%s}. %s. %s ' % (
                    entry.authors, entry.year, title, location, term)
                f.write(line + '\n\n')


def write_talks(talks_file='talks.csv'):
    if os.path.exists(talks_file):
        talks = pd.read_csv(talks_file, index_col=0)
        years = list(talks.year.unique())
        years.sort()
        years = years[::-1]

        term = '\\vspace{2mm}'
        lines = []
        for y in years:
            talks_year = talks.query('year == %s' % y)
            lines.append('%s: %s' % (y, ','.join(list(talks_year.place))))
        print('writing talks to talks.tex')
        with open('talks.tex', 'w') as f:
            f.write('\\section*{Invited addresses and colloquia (* - talks given virtually)}\n\\noindent\n\n')
            for l in lines:
                f.write(l + term + '\n\n')


def get_doi_from_pmid(pmid, pubs):
    # utility function
    for p in pubs:
        if pubs[p]['pmid'] == pmid:
            print('%s: %s' % (pubs[p]['pmid'], p))
