"""
class for a researcher
"""

import os
import json
import requests
import scholarly
import pypatent
from .orcid import get_dois_from_orcid_record
from .pubmed import get_pubmed_data
from .publication import JournalArticle, get_random_hash
from .crossref import get_crossref_records, parse_crossref_record
from .latex import render_latex


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
        self.publications = None
        self.rendered_latex = None

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
        self.orcid_dois = get_dois_from_orcid_record(self.orcid_data)

    def get_pubmed_data(self):
        self.pubmed_data = get_pubmed_data(self.query, self.email)
        print('retrieved %d full pubmed records' % len(self.pubmed_data['PubmedArticle']))

    def get_google_scholar_record(self):
        search_query = scholarly.scholarly.search_author(
            ' '.join([self.firstname, self.lastname]))
        self.gscholar_data = next(search_query).fill()

    def make_publication_records(self):
        # test pubmed
        self.get_pubmed_data()
        pubmed_dois = []
        self.publications = {}
        for r in self.pubmed_data['PubmedArticle']:
            pub = JournalArticle()
            pub.from_pubmed(r)
            pub.format_reference_latex()
            pub.hash = pub.get_pub_hash()
            self.publications[pub.DOI] = pub
            # keep track of pubmed DOIs so that we
            # don't overwrite with crossref
            pubmed_dois.append(pub.DOI)

        if self.orcid_data is None:
            self.get_orcid_data()
        if self.orcid_dois is None:
            self.get_orcid_dois()
        print('found %d  ORCID dois' % len(self.orcid_dois))

        # load orcid pubs using crossref
        self.crossref_data = get_crossref_records(self.orcid_dois)
        print('found %d crossref records' % len(self.crossref_data))

        for c in self.crossref_data:
            d = parse_crossref_record(self.crossref_data[c])
            if d is not None:
                p = JournalArticle()
                p.from_dict(d)
                # p.format_reference_latex()
                p.hash = p.get_pub_hash()
                if p.DOI not in pubmed_dois:
                    self.publications[p.DOI] = p
        print('found %d additional pubs from ORCID via crossref' % (len(self.publications) - len(pubmed_dois)))

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

            self.serialized['publications'] = {}
            for k in self.publications:
                if self.publications[k].hash is None:
                    self.publications[k].hash = get_random_hash()
                self.serialized['publications'][self.publications[k].hash] = self.publications[k].to_json()

    def to_json(self, filename):
        if self.serialized is None:
            self.serialize()
        with open(filename, 'w') as f:
            json.dump(self.serialized, f)

    def render_latex_cv(self):
        self.rendered_latex = render_latex(self)
