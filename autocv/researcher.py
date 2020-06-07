"""
class for a researcher
"""

import os
import json
from pprint import pprint
import requests
from Bio import Entrez
import scholarly


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


class Researcher:

    def __init__(self, params_file='params.json'):
        self.load_params(params_file)
        self.orcid_data = None
        self.orcid_dois = None
        self.pubmed_data = None
        self.crossref_data = None
        self.gscholar_data = None
        self.patent_data = None


    def load_params(self, params_file):
        if os.path.exists(params_file):
            with open(params_file) as f:
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

    
if __name__ == '__main__':
    r = Researcher('../tests/params.json')
    pprint(vars(r))
    r.get_orcid_data()
    r.get_pubmed_data()
