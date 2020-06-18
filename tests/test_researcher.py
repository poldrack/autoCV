"""
tests for researcher class
using params for Russ Poldrack as test case
"""

import pytest
import os
import pkg_resources
from autocv.researcher import Researcher
from autocv.orcid import get_orcid_education
from autocv.orcid import get_orcid_funding
from autocv.orcid import get_orcid_employment
from autocv.orcid import get_orcid_memberships
from autocv.orcid import get_orcid_service


# test fixtures to share across tests
@pytest.fixture(scope="session")
def researcher():
    param_file = pkg_resources.resource_filename('autocv', 'testdata/params.json')
    r = Researcher(param_file)
    return(r)


@pytest.fixture(scope="session")
def jsonfile(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("test.json")
    return(fn)


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
    assert researcher.gscholar_data.hindex >= 112


def test_researcher_get_patents(researcher):
    researcher.get_patents()
    # use number as of June 7, 2020
    assert len(researcher.patent_data) >= 1


def test_orcid_get_education(researcher):
    education = get_orcid_education(researcher.orcid_data)
    assert education.shape[0] >= 3


def test_orcid_get_funding(researcher):
    funding = get_orcid_funding(researcher.orcid_data)
    assert funding.shape[0] >= 23


def test_orcid_get_employment(researcher):
    employment = get_orcid_employment(researcher.orcid_data)
    assert employment.shape[0] >= 12


def test_orcid_get_memberships(researcher):
    memberships = get_orcid_memberships(researcher.orcid_data)
    assert memberships.shape[0] >= 4


def test_orcid_get_service(researcher):
    service = get_orcid_service(researcher.orcid_data)
    assert service.shape[0] >= 32


def test_make_publication_records(researcher):
    researcher.make_publication_records()
    assert len(researcher.publications) >= 280


def test_serialize(researcher):
    researcher.serialize()
    assert isinstance(researcher.serialized, dict)


def test_to_json(researcher, jsonfile):
    researcher.to_json(jsonfile)
    assert os.path.exists(jsonfile)


def test_from_json(researcher, jsonfile):
    researcher_tmp = Researcher(researcher.param_file)
    researcher_tmp.from_json(jsonfile)
    # spot check
    assert len(researcher.pubmed_data['PubmedArticle']) ==\
        len(researcher_tmp.pubmed_data['PubmedArticle'])
