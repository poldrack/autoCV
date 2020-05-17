"""
generate a LaTeX CV using PubMed and ORCID along with other resources

Russ Poldrack, May 2020
"""

import pandas as pd
from Bio import Entrez
import json
import os
from pprint import pprint
import random
import string
import requests

## settings

drop_corrigenda = True

# load parameters

def get_params(params_file='params.json'):
    if os.path.exists(params_file):
        with open(params_file) as f:
            params = json.load(f)
    else:
        raise FileNotFoundError('Please create a json file called params.json containing the fields email (with your email address), orcid (with your ORCID id) and query (with your pubmed query)- see documentation for help')
    assert "orcid" in params
    assert "email" in params
    assert "query"
    return(params)

def get_orcid_data(id):
    resp = requests.get("http://pub.orcid.org/%s"%id,
                    headers={'Accept':'application/orcid+json'})
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
        ctr +=1

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
        ctr +=1
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
        ctr +=1

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



def get_orcid_pubs(orcid_data):
    pubs = {}
    for g in orcid_data['activities-summary']['works']['group']:
        for p in g['work-summary']:
            if not p['type'] == 'JOURNAL_ARTICLE':
                continue
            title=p['title']['title']['value']
            doi=None
            for eid in p['external-ids']['external-id']:
                if eid['external-id-type']== 'doi':
                    doi = eid['external-id-value']
            if doi is not None:
                pubs[doi]=title
    return(pubs)


def render_education(education_df):
    with open('education.tex', 'w') as f:
        for i in education_df.index:

            f.write('\\textit{%s-%s}: %s (%s), %s, %s\n\n' %(
                education_df.loc[i, 'start_date'],
                education_df.loc[i, 'end_date'],
                education_df.loc[i, 'degree'],
                education_df.loc[i, 'dept'],
                education_df.loc[i, 'institution'],
                education_df.loc[i, 'city'],
            ))

def render_employment(employment_df):
    with open('employment.tex', 'w') as f:
        for i in employment_df.index:
            if employment_df.loc[i, 'dept'] is None:
                dept = ''
            else:
                dept = ' (%s)' % employment_df.loc[i, 'dept']

            f.write('\\textit{%s-%s}: %s%s, %s\n\n' %(
                employment_df.loc[i, 'start_date'],
                employment_df.loc[i, 'end_date'],
                employment_df.loc[i, 'role'],
                dept,
                employment_df.loc[i, 'institution'],
            ))

def render_distinctions(distinctions_df):
    with open('distinctions.tex', 'w') as f:
        for i in distinctions_df.index:
            f.write('\\textit{%s}: %s, %s\n\n' %(
                distinctions_df.loc[i, 'start_date'],
                distinctions_df.loc[i, 'title'],
                distinctions_df.loc[i, 'organization'],
            ))

def render_memberships(memberships_df):
     with open('memberships.tex', 'w') as f:
        memberships = ', '.join(memberships_df.organization)
        f.write('%s\n\n' % memberships)
   
# read publications from PubMed
def get_pubmed_pubs(params):
    Entrez.email = params['email']
    print(f'using {Entrez.email} for Entrez service')
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
    print('found %d full records' % len(records['PubmedArticle']))
    return(records)

#if __name__ == "__main__":

params = get_params()

# get orcid data
print('reading data for ORCID id:', params['orcid'])
orcid_data = get_orcid_data(params['orcid'])
education_df = get_orcid_education(orcid_data)
employment_df = get_orcid_employment(orcid_data)
distinctions_df = get_orcid_distinctions(orcid_data)
memberships_df = get_orcid_memberships(orcid_data)

pubmed_records = get_pubmed_pubs(params)

# load various links
## OSF repositories
osf_links = {}
if os.path.exists('osf_links.csv'):
    osf_df = pd.read_csv('osf_links.csv', index_col=0)
    for i in osf_df.index:
        osf_links[i] = osf_df.loc[i, 'url']
    print('OSF links:')
    pprint(osf_links)
else:
    print('OSF links file (osf_links.csv) not found, skipping...')

## Code links
code_links = {}
if os.path.exists('code_links.csv'):
    code_df = pd.read_csv('code_links.csv', index_col=0)
    for i in code_df.index:
        code_links[i] = code_df.loc[i, 'url']
    print('Code links:')
    pprint(code_links)

## Data links
data_links = {}
if os.path.exists('data_links.csv'):
    data_df = pd.read_csv('data_links.csv', index_col=0)
    for i in data_df.index:
        data_links[i] = data_df.loc[i, 'url']
    print('Data links:')
    pprint(data_links)

#
## Parse pubmed records to get reference info

pubs = {}

for i in records['PubmedArticle']:
    pmid = int(i['MedlineCitation']['PMID'])
    if 'Year' in i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
        year = i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
    elif 'MedlineDate' in i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
        year = i['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'].split(' ')[0]
    else:
        print('no year - skipping', pmid)

    if year not in pubs:
        pubs[year] = {}
    pubs[year][pmid] = {}

    authorlist = []
    if 'AuthorList' in i['MedlineCitation']['Article']:
        for author in i['MedlineCitation']['Article']['AuthorList']:
            if 'LastName' in author and 'Initials' in author:
                authorlist.append(' '.join([author['LastName'], author['Initials']]))
        pubs[year][pmid]['authors'] = ', '.join(authorlist)
    if 'Volume' in i['MedlineCitation']['Article']['Journal']['JournalIssue']:
        pubs[year][pmid]['volume'] = i['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
    else:
        pubs[year][pmid]['volume'] = ''

    pubs[year][pmid]['title'] = i['MedlineCitation']['Article']['ArticleTitle']

    if 'Pagination' in i['MedlineCitation']['Article']:
        pubs[year][pmid]['pages'] = i['MedlineCitation']['Article']['Pagination']['MedlinePgn']

    pubs[year][pmid]['journal'] = i['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
    # try to get PMC
    for j in i['PubmedData']['ArticleIdList']:
        if j.attributes['IdType'] == 'pmc':
            pubs[year][pmid]['PMC'] = str(j)
        if j.attributes['IdType'] == 'doi':
            pubs[year][pmid]['DOI'] = str(j)

    if pmid in osf_links:
        pubs[year][pmid]['OSF'] = osf_links[pmid]

    if pmid in code_links:
        pubs[year][pmid]['code'] = code_links[pmid]

    if pmid in data_links:
        pubs[year][pmid]['data'] = data_links[pmid]


# load addtional publications from csv
if os.path.exists('additional_pubs.csv'):
    addpubs = pd.read_csv('additional_pubs.csv')
    addpubs = addpubs.fillna('')
    for i in addpubs.index:
        year = str(addpubs.loc[i, 'year'])
        # make a random string to stand in for pmid
        pmid = ''.join(random.choice(
            string.ascii_lowercase) for i in range(8))
        pubs[year][pmid] = {}
        for c in addpubs.columns:
            pubs[year][pmid][c] = addpubs.loc[i, c]

# load chapters
if os.path.exists('chapters.csv'):
    chapters = pd.read_csv('chapters.csv')
    chapters = chapters.fillna('')
    for i in chapters.index:
        year = str(chapters.loc[i, 'year'])
        # make a random string to stand in for pmid
        pmid = ''.join(random.choice(
            string.ascii_lowercase) for i in range(8))
        pubs[year][pmid] = {}
        pubs[year][pmid] = {'isChapter': True}
        for c in chapters.columns:
            pubs[year][pmid][c] = chapters.loc[i, c]

# load books
if os.path.exists('books.csv'):
    books = pd.read_csv('books.csv')
    books = books.fillna('')
    for i in books.index:
        year = str(books.loc[i, 'year'])
        # make a random string to stand in for pmid
        pmid = ''.join(random.choice(
            string.ascii_lowercase) for i in range(8))
        pubs[year][pmid] = {'isBook': True}
        for c in books.columns:
            pubs[year][pmid][c] = books.loc[i, c]


## convert to latex

# get list of years
years = list(pubs.keys())
years.sort(reverse=True)

latex_lines = []
for year in years:
    latex_lines.append('\\subsection*{%s}' % year)
    # first get alphabetical order of PMIDs
    pmid_df = pd.DataFrame({'author': ''}, index=list(pubs[year].keys()))
    for pmid in pubs[year]:
        pmid_df.loc[pmid, 'author'] = pubs[year][pmid]['authors']
    pmid_df.sort_values('author', inplace=True)
    pmids_sorted = list(pmid_df.index)

   

    # get the reference line for each
    for pmid in pmids_sorted:
        if pubs[year][pmid]['title'].find('Corrigend') > -1:
            continue
        if len(pubs[year][pmid]['authors'].split(',')) > 10:
            authorlist = pubs[year][pmid]['authors'].split(',')[0] + ', et al.'
        else:
            authorlist = pubs[year][pmid]['authors']
        if 'isBook' in pubs[year][pmid]:  # treat books separately
            line = authorlist +\
                ' (%s). ' % year +\
                ' \\textit{%s} ' % pubs[year][pmid]['title'] +\
                pubs[year][pmid]['publisher'].strip(' ')
            line += '.'
        elif 'isChapter' in pubs[year][pmid]:  # treat chapters separately
            line = authorlist +\
                ' (%s). ' % year +\
                pubs[year][pmid]['title'] +\
                'In %s (Eds.),  \\textit{%s.} %s.' % (
                    pubs[year][pmid]['editors'],
                    pubs[year][pmid]['booktitle'],
                    pubs[year][pmid]['publisher'].strip(' '))
        else:
            line = authorlist +\
                ' (%s). ' % year +\
                pubs[year][pmid]['title'] +\
                ' \\textit{%s' % pubs[year][pmid]['journal']
            if pubs[year][pmid]['volume'] != '':
                line += ', %s}' % pubs[year][pmid]['volume']
            else:
                line += '}'
            if 'pages' in pubs[year][pmid]:
                if len(pubs[year][pmid]['pages']) > 0:
                    line += ', %s' % pubs[year][pmid]['pages']
            line += '.'

        if 'PMC' in pubs[year][pmid]:
            line += ' \\href{https://www.ncbi.nlm.nih.gov/pmc/articles/%s}{OA}' % pubs[year][pmid]['PMC']

        if 'OSF' in pubs[year][pmid]:
            line += ' \\href{%s}{OSF}' % pubs[year][pmid]['OSF']

        if 'code' in pubs[year][pmid]:
            line += ' \\href{%s}{Code}' % pubs[year][pmid]['code']

        if 'data' in pubs[year][pmid]:
            line += ' \\href{%s}{Data}' % pubs[year][pmid]['data']

        if 'DOI' in pubs[year][pmid]:
            if len(pubs[year][pmid]['DOI']) > 0:
                line += ' \\href{http://dx.doi.org/%s}{DOI}' % pubs[year][pmid]['DOI']

        line += ' \\vspace{2mm}'
        latex_lines.append(line)

# write pubs to file
print('writing pubs to pubs.tex')
with open('pubs.tex', 'w') as f:
    for l in latex_lines:
        f.write(l.replace('_', '\_') + '\n\n') # noqa


# ### load presentations and write to latex

# +
if os.path.exists('presentations.csv'):
    presentations = pd.read_csv('presentations.csv', index_col=0)
    presentations = presentations.sort_values('year', ascending=False)

    print('writing presentations to presentations.tex')
    with open('presentations.tex', 'w') as f:
        for i in presentations.index:
            entry = presentations.loc[i, :]
            title = entry.title.strip('.')
            location = entry.location.strip(' ').strip('.')
            term = '\\vspace{2mm}'
            line = '%s (%s). \\emph{%s}. %s. %s ' % (
                entry.authors, entry.year, title, location, term)
            f.write(line + '\n\n')


# separately print colloquium talks
if os.path.exists('talks.csv'):
    talks = pd.read_csv('talks.csv', index_col=0)
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
        for l in lines:
            f.write(l + term + '\n\n')
