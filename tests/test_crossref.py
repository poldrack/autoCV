"""
test crossref utilty functions in crossref.py
"""

import pytest
from autocv.crossref import get_crossref_records, parse_crossref_record


@pytest.fixture(scope="session")
def dois():
    return(['10.7554/elife.53498', '10.1016/j.neuron.2020.02.019'])


@pytest.fixture(scope="session")
def crossref_records(dois):
    crossref_records = get_crossref_records(dois)
    return(crossref_records)


def test_get_crossref_records(crossref_records, dois):
    assert len(crossref_records) == len(dois)
    for d in dois:
        assert d in crossref_records

def test_parse_crossref_records(crossref_records):
    for record in crossref_records:
        parsed_record = parse_crossref_record(crossref_records[record])