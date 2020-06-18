"""
tests for researcher class
using params for Russ Poldrack as test case
"""

import pytest
import pkg_resources
from autocv.researcher import Researcher
from autocv.latex import LatexCV


# test fixtures to share across tests
@pytest.fixture(scope="session")
def latexcv():
    param_file = pkg_resources.resource_filename('autocv', 'testdata/params.json')
    r = Researcher(param_file)
    r.get_orcid_data()
    r.get_orcid_dois()
    r.get_pubmed_data()
    r.get_google_scholar_record()
    r.get_patents()
    r.make_publication_records()
    return(LatexCV(r))


#  really just smoke tests
def test_latex_cv_class(latexcv):
    assert latexcv is not None


def test_latex_cv_load_template_files(latexcv):
    latexcv.load_template_files()
    assert latexcv.front is not None
    assert latexcv.back is not None


def test_latex_cv_render_latex(latexcv):
    latexcv.render_latex()
    for section in latexcv.sections_to_write:
        assert hasattr(latexcv, section)


def test_latex_cv_write_latex(latexcv, tmpdir_factory):
    outfile = tmpdir_factory.mktemp("data").join("test.tex")
    latexcv.write_latex(outfile)
