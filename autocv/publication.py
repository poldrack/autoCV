"""
class for publications
"""

import random
import string
import hashlib
import json

from .crossref import get_crossref_records, parse_crossref_record
from .pubmed import parse_pubmed_record


def get_random_hash(length=16):
    return(''.join(random.choice(string.ascii_lowercase) for i in range(length)))


def serialize_pubs_to_json(pubs, outfile):
    """
    save a list of publications to json

    parameters:
    -----------
    pubs: a list of Publication objects
    outfile: string, filename to save to
    """

    # first combine into a single dictionary
    pubdict = {}
    for p in pubs:
        if p.hash in pubdict:
            print('WARNING: hash collision')
            p.hash = p.hash + get_random_hash(4)
        pubdict[p.hash] = vars(p)
    with open(outfile, 'w') as f:
        json.dump(pubdict, f)
    return(pubdict)


def load_pubs_from_json(infile):
    pubdict = {}
    with open(infile) as f:
        pubdict = json.load(f)
    return(pubdict)


class Publication:
    """
    """

    type = 'generic'

    def __init__(self, title=None, year=None, authors=None):

        # set up general feature attributes
        self.title = title
        self.year = year
        self.authors = authors
        self.hash = None

    def get_pub_hash(self, digest_size=8):
        """
        create a hash from the title, year, and authors
        - used for finding duplicates
        """
        if self.title is None:
            print('reference must first be loaded')
        else:
            pubstr = '-'.join([str(i) for i in [self.title, self.year, self.authors]])
            self.hash = hashlib.blake2b(pubstr.lower().encode('utf-8'), digest_size=digest_size).hexdigest()
    
    def to_json(self):
        return(vars(self))


class JournalArticle(Publication):

    type = 'journal-article'

    def __init__(self, title=None, year=None, authors=None,
                 journal=None, volume=None, page=None, DOI=None):
        super().__init__(title, year, authors)

        self.journal = journal
        self.volume = volume
        self.page = page
        self.DOI = DOI
        self.PMC = None
        self.PMID = None
        self.links = {}
        self.reference = None
        self.source = None
        self.type = 'journal-article'

    def format_reference_latex(self):
        if self.title is None:
            print('reference must be loaded before formatting')
            return

        line = self.authors +\
            ' (%d). ' % self.year +\
            self.title +\
            ' \\textit{%s' % self.journal

        line += ', %s}' % self.volume if self.volume is not None else '}'
        if self.page is not None and len(self.page) > 0:
            line += ', %s' % self.page
        line += '.'
        self.reference = line

    def from_dict(self, pubdict):
        for k in pubdict:
            if hasattr(self, k):
                setattr(self, k, pubdict[k])

    def from_pubmed(self, pubmed_record):
        parsed_record = parse_pubmed_record(pubmed_record)
        self.source = 'Pubmed'
        for k in parsed_record:
            setattr(self, k, parsed_record[k])


if __name__ == "__main__":
    rsrchr = Researcher('../tests/params.json')

    # test pubmed
    pubmed_records = rsrchr.get_pubmed_records('poldrack-r', 'poldrack@stanford.edu')
    pubmed_publications = []
    pubmed_dois = []
    for r in pubmed_records['PubmedArticle']:
        pub = JournalArticle()
        pub.from_pubmed(r)
        pub.format_reference_latex()
        pub.hash = pub.get_pub_hash()
        pubmed_publications.append(pub)
        pubmed_dois.append(pub.DOI)

    # test orcid
    orcid_data = rsrchr.get_orcid_data()
    orcid_dois = rsrchr.get_orcid_dois()
    print('found %d  ORCID dois' % len(orcid_dois))

    # load orcid pubs using crossref
    crossref_records = get_crossref_records(orcid_dois)
    print('found %d crossref records' % len(crossref_records))

    crossref_pubs = []
    for c in crossref_records:
        d = parse_crossref_record(crossref_records[c])
        if d is not None:
            p = JournalArticle()
            p.from_dict(d)
            # p.format_reference_latex()
            p.hash = p.get_pub_hash()
            if p.DOI not in pubmed_dois:
                crossref_pubs.append(p)
    print('found %d additional pubs from ORCID via crossref' % len(crossref_pubs))

    # test saving
    pubs = pubmed_publications + crossref_pubs
    pubs_dict = serialize_pubs_to_json(pubs, 'test.json')
    pubs_retrieved = load_pubs_from_json('test.json')

    for i in range(len(pubs)):
        assert pubs[i].__dict__ == pubs_retrieved[i].__dict__
