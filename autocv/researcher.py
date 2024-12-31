"""
class for a researcher
"""

import os
import json
import requests
import scholarly
from scholarly import MaxTriesExceededException
import pypatent
from urllib.error import HTTPError

from .orcid import get_dois_from_orcid_record
from .pubmed import get_pubmed_data
from .publication import JournalArticle, Book, BookChapter
from .crossref import get_crossref_records, parse_crossref_record
from .utils import get_additional_pubs_from_csv, CustomJSONEncoder, get_random_hash, drop_excluded_pubs


class Researcher:

    def __init__(self, param_file='params.json', basedir=None, format='latex'):
        self.param_file = param_file
        self.load_params(param_file)
        self.basedir = os.path.dirname(param_file) if basedir is None else basedir
        self.format = format
        fields_to_add = [
            'orcid_data',
            'orcid_dois',
            'pubmed_data',
            'crossref_data',
            'gscholar_data',
            'patent_data',
            'serialized',
            'publications',
        ]
        for field in fields_to_add:
            if not hasattr(self, field):
                setattr(self, field, None)

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

    def get_orcid_data(self, timeout=60):
        orcid_url = "https://pub.orcid.org/v3.0/%s" % self.orcid
        print('using ORCID URL:', orcid_url)
        resp = requests.get(orcid_url,
                            headers={'Accept': 'application/vnd.orcid+json'},
                            timeout=timeout)
        self.orcid_data = resp.json()
        if 'error-code' in self.orcid_data:
            raise ValueError(
                f"problem accessing ORCID: {self.orcid_data['developer-message']}")

    def get_orcid_dois(self):
        if self.orcid_data is None:
            self.get_orcid_data()
        self.orcid_dois = get_dois_from_orcid_record(self.orcid_data)

    def get_pubmed_data(self):
        self.pubmed_data = get_pubmed_data(self.query, self.email)
        print('retrieved %d full pubmed records' % len(self.pubmed_data['PubmedArticle']))

    def get_google_scholar_record(self):
        try:
            search_query = scholarly.scholarly.search_author(
                ' '.join([self.firstname, self.lastname]))
            query_resp = next(search_query)
            self.gscholar_data = scholarly.scholarly.fill(query_resp)
        except MaxTriesExceededException:
            print('problem accessing google scholar')


    def make_publication_records(self, use_exclusions=True):
        # test pubmed
        self.get_pubmed_data()
        pubmed_dois = []
        self.publications = {}
        for r in self.pubmed_data['PubmedArticle']:
            pub = JournalArticle()
            pub.from_pubmed(r)

            pub.format_reference(format=self.format)
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
                # skip existing pubmed records and preprints
                if d['DOI'] in pubmed_dois:
                    continue

                if d['type'] in ['journal-article', 'proceedings-article']:
                    p = JournalArticle()
                elif d['type'] in ['book', 'monograph']:
                    p = Book()
                elif d['type'] == 'book-chapter':
                    p = BookChapter()
                else:
                    continue

                p.from_dict(d)
                if hasattr(p, 'DOI'):
                    id = p.DOI
                elif hasattr(p, 'ISBN'):
                    id = p.ISBN
                else:
                    id = get_random_hash()

                self.publications[id] = p
        if use_exclusions:
            self.publications = drop_excluded_pubs(self.publications)

        print('found %d additional pubs from ORCID via crossref' % (len(self.publications) - len(pubmed_dois)))

        additional_pubs_file = os.path.join(
            self.basedir, 'additional_pubs.csv'
        )
        additional_pubs = get_additional_pubs_from_csv(additional_pubs_file)
        for pub in additional_pubs:
            if additional_pubs[pub]['type'] in ['journal-article', 'proceedings-article']:
                self.publications[pub] = JournalArticle()
            elif additional_pubs[pub]['type'] in ['book', 'monograph']:
                self.publications[pub] = Book()
            elif additional_pubs[pub]['type'] == 'book-chapter':
                self.publications[pub] = BookChapter()
            else:
                print('skipping unknown type', additional_pubs[pub]['type'])
                continue
            self.publications[pub].from_dict(additional_pubs[pub])

    def get_patents(self):
        try:
            # this will fail if there are no patents
            results = pypatent.Search(self.lastname).as_list()
        except AttributeError:
            print('pypatent failed to retrieve patents')
            results = []

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
                # print('ingesting', k)
                if k == 'publications':
                    self.publications = {}
                    for pub in serialized[k]:
                        if serialized[k][pub]['type'] in ['journal-article', 'proceedings-article']:
                            self.publications[pub] = JournalArticle()
                        elif serialized[k][pub]['type'] in ['book', 'monograph']:
                            self.publications[pub] = Book()
                        elif serialized[k][pub]['type'] == 'book-chapter':
                            self.publications[pub] = BookChapter()
                        else:
                            print('skipping unknown type', serialized[k][pub]['type'])
                            continue
                        self.publications[pub].from_dict(serialized[k][pub])
                else:
                    setattr(self, k, serialized[k])

    def serialize_publications(self):
        self.serialized['publications'] = {}
        for k, pubinfo_orig in self.publications.items():
            pubinfo = pubinfo_orig.to_json()
            if len(pubinfo) == 0:
                print('skipping', k)
                continue
            else:
                self.serialized['publications'][k] = pubinfo

    def serialize(self):
        self.serialized = {}
        self_dict = self.__dict__.copy()
        if 'gscholar_data' in self_dict and self_dict['gscholar_data'] is not None and 'hindex' in self_dict['gscholar_data']:
            self.serialized['gscholar_data'] = {
                'hindex': self_dict['gscholar_data']['hindex']}
        
        self.serialize_publications()

    def to_json(self, filename):
        if self.serialized is None:
            self.serialize()
        with open(filename, 'w') as f:
            json.dump(self.serialized, f, cls=CustomJSONEncoder,
                      indent=4)
