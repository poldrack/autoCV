"""
class for publications
"""

import hashlib
import json

from .utils import get_random_hash
from .pubmed import parse_pubmed_record


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


def shorten_authorlist(authors, maxlen=10, n_to_show=3):
    authors_split = authors.split(',')
    if len(authors_split) > maxlen:
        authors = ','.join(authors_split[:n_to_show]) + ' et al.'
    return authors


def load_pubs_from_json(infile):
    pubdict = {}
    with open(infile) as f:
        pubdict = json.load(f)
    return(pubdict)


class Publication:
    """
    """

    type = 'generic'

    def __init__(self, title=None, year=None, authors=None, etalthresh=10):

        # set up general feature attributes
        self.title = title
        self.year = year
        self.authors = authors
        self.etalthresh = etalthresh
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

    def from_dict(self, pubdict):
        for k in pubdict:
            if hasattr(self, k):
                setattr(self, k, pubdict[k])

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

    def format_reference_latex(self, etalthresh=10, etalnum=3):

        if self.title is None:
            print('reference must be loaded before formatting')
            return
        authors_shortened = shorten_authorlist(self.authors, etalthresh, etalnum)

        line = authors_shortened +\
            ' (%d). ' % self.year +\
            self.title +\
            ' \\textit{%s' % self.journal

        line += ', %s}' % self.volume if self.volume is not None else '}'
        if self.page is not None and len(self.page) > 0:
            line += ', %s' % self.page
        line += '.'
        return(line)

    def from_pubmed(self, pubmed_record):
        parsed_record = parse_pubmed_record(pubmed_record)
        self.source = 'Pubmed'
        for k in parsed_record:
            setattr(self, k, parsed_record[k])


class BookChapter(Publication):

    type = 'book-chapter'

    def __init__(self, title=None, year=None, authors=None,
                 journal=None, page=None, ISBN=None,
                 publisher=None, editors=None):
        super().__init__(title, year, authors)

        self.journal = journal
        self.page = page
        self.ISBN = ISBN
        self.links = {}
        self.reference = None
        self.source = None
        self.publisher = publisher
        self.editors = editors

    def format_reference_latex(self, etalthresh=None, etalnum=None):
        if self.title is None:
            print('reference must be loaded before formatting')
            return

        page_string = ''
        if hasattr(self, 'page') and self.page is not None and len(self.page) > 0:
            page_string = '(p. %s). ' % self.page
        return self.authors +\
            ' (%s). ' % self.year +\
            self.title.strip('.') +\
            '. In \\textit{%s.} %s%s.' % (
                self.journal,
                page_string,
                self.publisher.strip(' '))


class Book(Publication):

    type = 'book'

    def __init__(self, title=None, year=None, authors=None,
                 page=None, ISBN=None,
                 publisher=None, editors=None):
        super().__init__(title, year, authors)

        self.page = page
        self.ISBN = ISBN
        self.links = {}
        self.reference = None
        self.source = None
        self.publisher = publisher
        self.editors = editors

    def format_reference_latex(self, etalthresh=None, etalnum=None):
        if self.title is None:
            print('reference must be loaded before formatting')
            return
        line = self.authors +\
            ' (%s). ' % self.year +\
            ' \\textit{%s}. ' % self.title.strip(' ').strip('.') + \
            self.publisher.strip(' ')
        line += '.'
        return(line)
