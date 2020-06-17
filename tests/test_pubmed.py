"""
tests for pubmed utility functions in pubmed.py
"""

import pytest
from Bio import Entrez
from autocv.pubmed import parse_pubmed_pubs, get_pubmed_doi


@pytest.fixture(scope="session")
def pubmed_data():
    pmids = [32425159, 32272064]
    Entrez.email = 'testing@autocv.org'
    handle = Entrez.efetch(
        db="pubmed",
        id=",".join(['%d' % i for i in pmids]),
        retmax=1000, retmode="xml")
    return(Entrez.read(handle))


def test_parse_pubmed_pubs(pubmed_data):
    pubs = parse_pubmed_pubs(pubmed_data)
    assert isinstance(pubs, dict)
    assert len(pubs) == 2
    for i, pub in enumerate(pubmed_data['PubmedArticle']):
        doi = get_pubmed_doi(pubmed_data['PubmedArticle'][i])
        assert doi in pubs


