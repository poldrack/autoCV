"""
functions to access various APIs
"""

import requests
from Bio import Entrez
from crossref.restful import Works
import pypatent
import scholarly


def get_orcid_data(id):
    resp = requests.get("http://pub.orcid.org/%s" % id,
                        headers={'Accept': 'application/orcid+json'})
    orcid_data = resp.json()
    return(orcid_data)


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


def get_pubmed_records(query, email=None):
    if email is not None:
        Entrez.email = email
        print(f'using {email} for Entrez service')
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
    print('retrieved %d full pubmed records' % len(records['PubmedArticle']))
    return(records)


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


def get_google_scholar_record(firstname, lastname):
    search_query = scholarly.scholarly.search_author(' '.join([firstname, lastname]))
    author = next(search_query).fill()
    return(author)


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
