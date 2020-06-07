"""
tests related to Publication class
"""

import pytest
from Bio import Entrez
import pkg_resources
from autocv.researcher import Researcher
from autocv.publication import Publication, JournalArticle
from autocv.crossref import get_crossref_records, parse_crossref_record

# test fixtures to share across tests
@pytest.fixture(scope="session")
def publication():
    p = Publication()
    return(p)


@pytest.fixture(scope="session")
def pubmed_data():
    pmids = [32425159, 32272064]
    Entrez.email = 'testing@autocv.org'
    handle = Entrez.efetch(
        db="pubmed",
        id=",".join(['%d' % i for i in pmids]),
        retmax=1000, retmode="xml")
    return(Entrez.read(handle))


def test_publication_class():
    pub = Publication()
    assert pub is not None


def test_article_class():
    pub = JournalArticle()
    assert pub is not None


def test_pubmed_data(pubmed_data):
    assert pubmed_data is not None


def test_article_from_pubmed(pubmed_data):
    for record in pubmed_data['PubmedArticle']:
        j = JournalArticle()
        j.from_pubmed(record)
        assert j.title == record['MedlineCitation']['Article']['ArticleTitle']
