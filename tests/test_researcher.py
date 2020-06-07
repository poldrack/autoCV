"""
tests for researcher class
using params for Russ Poldrack as test case
"""

import pytest
import pkg_resources

from autocv.researcher import Researcher
from autocv.publication import Publication


# test fixtures to share across tests
@pytest.fixture(scope="session")
def publication():
    p = Publication()
    return(p)


@pytest.fixture(scope="session")
def researcher():
    param_file = pkg_resources.resource_filename('autocv', 'testdata/params.json')
    r = Researcher(param_file)
    return(r)


# tests
def test_researcher_class(researcher):
    assert researcher is not None


def test_researcher_get_pubmed_data(researcher):
    researcher.get_pubmed_data()
    assert 'PubmedArticle' in researcher.pubmed_data
    # use number of pubs as of June 7, 2020 as test spec
    assert len(researcher.pubmed_data['PubmedArticle']) >= 261


def test_researcher_get_orcid_data(researcher):
    researcher.get_orcid_data()
    assert 'activities-summary' in researcher.orcid_data


def test_researcher_get_orcid_dois(researcher):
    researcher.get_orcid_dois()
    # use number of pubs as of June 7, 2020 as test spec
    assert len(researcher.orcid_dois) >= 227


def test_researcher_get_google_scholar_record(researcher):
    researcher.get_google_scholar_record()
    # use h index as of June 7, 2020
    assert len(researcher.orcid_dois) >= 227


def test_publication_class(publication):
    assert publication is not None
