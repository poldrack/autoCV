"""
test functions from csv.py
"""

import pkg_resources
from autocv.csv import get_links_from_csv
from autocv.csv import add_additional_pubs_from_csv
from autocv.csv import get_teaching_from_csv
from autocv.csv import get_funding_from_csv


def test_get_links_from_csv():
    link_file = pkg_resources.resource_filename(
        'autocv', 'testdata/links.csv')
    links = get_links_from_csv(link_file)
    # spot checks
    assert 'Code' in links
    assert len(links['Data']) >= 17


def test_get_additional_links_from_csv():
    pubs = {}
    pub_file = pkg_resources.resource_filename(
        'autocv', 'testdata/additional_pubs.csv')
    pubs = add_additional_pubs_from_csv(pubs, pubfile=pub_file)
    # length as of June 2020
    assert len(pubs) >= 22


def test_get_teaching_from_csv():
    teaching_file = pkg_resources.resource_filename(
        'autocv', 'testdata/teaching.csv')
    teaching = get_teaching_from_csv(teaching_file)
    # length as of June 2020
    assert len(teaching['Undergraduate']) >= 6


def test_get_funding_from_csv():
    funding_file = pkg_resources.resource_filename(
        'autocv', 'testdata/funding.csv')
    funding = get_funding_from_csv(funding_file)
    # length as of June 2020
    assert funding.shape[0] >= 37
