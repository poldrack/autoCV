"""
functions to access various APIs
"""

from crossref.restful import Works
import pypatent
import scholarly


def get_crossref_records(dois):
    works = Works()
    crossref_records = {}
    print('searching crossref for all DOIs, this might take a few minutes...')
    for doi in dois:
        r = works.doi(doi)
        if r is not None:
            crossref_records[doi] = r
        else:
            print('missing crossref record for', doi)
    return(crossref_records)


def get_google_scholar_record(firstname, lastname):
    search_query = scholarly.scholarly.search_author(' '.join([firstname, lastname]))
    author = next(search_query).fill()
    return(author)


def get_patents(lastname, firstname):
    results = pypatent.Search(lastname).as_list()
    mypatents = []
    for r in results:
        for i in r['inventors']:
            fn = i[0].split(' ')[0].lower()
            ln = i[1].lower()
            if fn == firstname.lower() and ln == lastname.lower():
                mypatents.append(r)
    return(mypatents)
