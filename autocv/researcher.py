"""
class for a researcher
"""

import os
import json
from pprint import pprint
import requests
from Bio import Entrez
import scholarly
import pypatent


def get_params(param_file='params.json'):
    if os.path.exists(param_file):
        with open(param_file) as f:
            params = json.load(f)
    else:
        raise FileNotFoundError('Please create a json file called params.json containing the fields email (with your email address), orcid (with your ORCID id) and query (with your pubmed query)- see documentation for help')
    required_fields = ['address', 'lastname', 'firstname', 'email', 'orcid', 'query', 'url', 'phone']
    for field in required_fields:
        assert field in params
    return(params)


class Researcher:

    def __init__(self, param_file='params.json'):
        self.param_file = param_file
        self.load_params(param_file)
        self.orcid_data = None
        self.orcid_dois = None
        self.pubmed_data = None
        self.crossref_data = None
        self.gscholar_data = None
        self.patent_data = None
        self.serialized = None

    def load_params(self, param_file):
        if os.path.exists(param_file):
            with open(param_file) as f:
                params = json.load(f)
        else:
            raise FileNotFoundError("""Please create a json file called params.json
                                       containing the fields email (with your email address), orcid (with your ORCID id)
                                       and query (with your pubmed query)- see documentation for help')
                                       """)
        for field in params:
            setattr(self, field, params[field])

    def get_orcid_data(self):
        resp = requests.get("http://pub.orcid.org/%s" % self.orcid,
                            headers={'Accept': 'application/orcid+json'})
        self.orcid_data = resp.json()

    def get_orcid_dois(self):
        if self.orcid_data is None:
            self.get_orcid_data()
        dois = []
        for g in self.orcid_data['activities-summary']['works']['group']:
            for p in g['work-summary']:
                doi = None
                for eid in p['external-ids']['external-id']:
                    if eid['external-id-type'] == 'doi':
                        doi = eid['external-id-value'].replace('http://dx.doi.org/', '')
                if doi is not None:
                    dois.append(doi.lower())
        self.orcid_dois = list(set(dois))

    def get_pubmed_data(self):
        Entrez.email = self.email
        print(f'using {self.email} for Entrez service')
        print('searching for', self.query)
        retmax = 1000
        handle = Entrez.esearch(db="pubmed", retmax=retmax, term=self.query)
        record = Entrez.read(handle)
        handle.close()
        pmids = [int(i) for i in record['IdList']]
        print('found %d matches' % len(pmids))

        # load full records
        handle = Entrez.efetch(db="pubmed", id=",".join(['%d' % i for i in pmids]),
                               retmax=retmax, retmode="xml")
        self.pubmed_data = Entrez.read(handle)
        print('retrieved %d full pubmed records' % len(self.pubmed_data['PubmedArticle']))

    def get_google_scholar_record(self):
        search_query = scholarly.scholarly.search_author(
            ' '.join([self.firstname, self.lastname]))
        self.gscholar_data = next(search_query).fill()

    def get_patents(self):
        results = pypatent.Search(self.lastname).as_list()
        self.patent_data = []
        for r in results:
            for i in r['inventors']:
                fn = i[0].split(' ')[0].lower()
                ln = i[1].lower()
                if fn == self.firstname.lower() and ln == self.lastname.lower():
                    self.patent_data.append(r)

    def from_json(self, filename):
        with open(filename, 'r') as f:
            serialized = json.load(f)
        for k in serialized.keys():
            if hasattr(self, k):
                setattr(self, k, serialized[k])

    def serialize(self):
        self.serialized = self.__dict__.copy()
        if 'gscholar_data' in self.serialized:
            # need to convert gscholar objects to dicts
            self.serialized['gscholar_data'] = self.serialized['gscholar_data'].__dict__

            coauthor_data = self.serialized['gscholar_data']['coauthors'].copy()
            self.serialized['gscholar_data']['coauthors'] = []
            for k in coauthor_data:
                self.serialized['gscholar_data']['coauthors'].append(k.__dict__)

            publication_data = self.serialized['gscholar_data']['publications'].copy()
            self.serialized['gscholar_data']['publications'] = []
            for k in publication_data:
                self.serialized['gscholar_data']['publications'].append(k.__dict__)

    def to_json(self, filename):
        if self.serialized is None:
            self.serialize()
        with open(filename, 'w') as f:
            json.dump(self.serialized, f)


if __name__ == '__main__':
    r = Researcher('autocv/testdata/params.json')
    pprint(vars(r))
    r.get_orcid_data()
    r.get_orcid_dois()
    r.get_pubmed_data()
    r.get_google_scholar_record()
    r.get_patents()
    testfile = 'test.json'
    r.to_json(testfile)
    r2 = Researcher(r.param_file)
    r2.from_json(testfile)
