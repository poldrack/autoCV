"""
class for publications
"""

import random
import string
import hashlib
import json

from .apis import get_crossref_records
from .researcher import Researcher


def parse_crossref_record(record, verbose=False, exclude_preprints=True,
                          exclude_books=True, exclude_translations=True):
    """
    extract fields from record
    do this here because these records span multiple publication types
    """
    pub = {'DOI': record['DOI']}
    if exclude_preprints and record['type'] == 'posted-content':
        if verbose:
            print('skipping preprint:', record['DOI'])
        return(None)
    # books seem to be goofy with crossref
    if exclude_books and record['type'] == 'book':
        if verbose:
            print('skipping book:', record['DOI'])
        return(None)
    if 'author' not in record:  # can happen for errata
        if verbose:
            print('skipping due to missing author:', record['DOI'])
        return(None)
    if 'translator' in record:
        if verbose:
            print('skipping translation:', record['DOI'])
        return(None)

    if isinstance(record['title'], list):
        record['title'] = record['title'][0]

    if record['title'].find('Corrigend') > -1:
        if verbose:
            print('skipping corrigendum:', record['DOI'])
        return(None)

    # don't replace pubmed info if it already exists
    for field in ['title', 'volume', 'page', 'type', 'publisher']:
        if field in record:
            f = record[field]
            if isinstance(f, list):
                f = f[0]
            pub[field] = f

    # filter out pages with n/a
    if 'page' in pub and pub['page'].find('n/a') > -1:
        del pub['page']

    # get the title
    if len(record['container-title']) > 0:
        pub['journal'] = record['container-title'][0]

    # date can show up in two different places!
    if 'published-print' in record:
        year = record['published-print']['date-parts'][0][0]
    elif 'journal-issue' in record:
        journal_issue = record['journal-issue']
        if 'published-print' in journal_issue:
            year = journal_issue['published-print']['date-parts'][0][0]
        else:
            year = journal_issue['published-online']['date-parts'][0][0]

        pub['year'] = int(year)

    # convert author list to pubmed format
    authors = []
    for author in record['author']:
        if 'given' not in author or 'family' not in author:
            continue
        given_split = author['given'].split(' ')
        if len(given_split) > 1:
            initials = ''.join([i[0] for i in given_split])
        else:
            initials = given_split[0][0]
        entry = '%s %s' % (author['family'], initials)
        authors.append(entry)
    pub['authors'] = ', '.join(authors)
    pub['source'] = 'Crossref'
    return(pub)


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
        self.source = 'Pubmed'
        self.PMID = int(pubmed_record['MedlineCitation']['PMID'])

        for j in pubmed_record['PubmedData']['ArticleIdList']:
            if j.attributes['IdType'] == 'doi':
                self.DOI = str(j).lower().replace('http://dx.doi.org/', '')
            if j.attributes['IdType'] == 'pmc':
                self.PMC = str(j)

        if self.DOI is None:
            print('no DOI found for', self.PMID)

        # get some other useful stuff while we are here
        # pubmed seems to have better info than crossref
        self.journal = pubmed_record['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        if 'Year' in pubmed_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            self.year = int(pubmed_record['MedlineCitation']['Article'][
                'Journal']['JournalIssue']['PubDate']['Year'])
        elif 'MedlineDate' in pubmed_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            self.year = int(pubmed_record['MedlineCitation']['Article'][
                'Journal']['JournalIssue']['PubDate']['MedlineDate'].split(' ')[0])
        if 'Volume' in pubmed_record['MedlineCitation']['Article']['Journal']['JournalIssue']:
            self.volume = pubmed_record['MedlineCitation']['Article'][
                'Journal']['JournalIssue']['Volume']
        self.title = pubmed_record['MedlineCitation']['Article']['ArticleTitle']

        if 'Pagination' in pubmed_record['MedlineCitation']['Article']:
            self.page = pubmed_record['MedlineCitation']['Article']['Pagination']['MedlinePgn']

        if 'AuthorList' in pubmed_record['MedlineCitation']['Article']:
            authorlist = [
                ' '.join([author['LastName'], author['Initials']])
                for author in pubmed_record['MedlineCitation']['Article'][
                    'AuthorList'
                ]
                if 'LastName' in author and 'Initials' in author
            ]

            self.authors = ', '.join(authorlist)


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
