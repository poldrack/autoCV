"""
functions to access crossref API
"""

from crossref.restful import Works


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


def parse_crossref_record(record, verbose=False, exclude_preprints=True,
                          exclude_books=True, exclude_translations=True):
    """
    extract fields from record
    do this here because these records span multiple publication types
    """
    if 'DOI' not in record:
        print('no DOI found in crossref record - skipping')
        print(record)
        return(None)
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
        elif 'published-online' in journal_issue:
            year = journal_issue['published-online']['date-parts'][0][0]
        else:
            print('problem getting year for:', pub)
            return(None)
    elif 'published-online' in record:
        year = record['published-online']['date-parts'][0][0]
    else:
        print('problem getting year for:', pub)
        return(None)

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


def process_crossref_records(crossref_records, pubs,
                             etal_thresh=10, exclude_preprints=True,
                             exclude_books=True, exclude_translations=True):
    # use DOI as dictionary keys
    for r in crossref_records:
        if exclude_preprints and crossref_records[r]['type'] == 'posted-content':
            continue
        # books seem to be goofy with crossref
        if exclude_books and crossref_records[r]['type'] == 'book':
            continue
        if 'author' not in crossref_records[r]:  # can happen for errata
            print('no author for ', r)
            continue
        if 'translator' in crossref_records[r]:
            continue

        if isinstance(crossref_records[r]['title'], list):
            crossref_records[r]['title'] = crossref_records[r]['title'][0]
        if crossref_records[r]['title'].find('Corrigend') > -1:
            continue

        if r not in pubs:
            pubs[r] = {}

        # don't replace pubmed info if it already exists
        for field in ['title', 'volume', 'page', 'type', 'publisher']:
            if field in crossref_records[r] and field not in pubs[r]:
                f = crossref_records[r][field]
                if isinstance(f, list):
                    f = f[0]
                pubs[r][field] = f

        # get the title
        if 'journal' not in pubs[r] and len(crossref_records[r]['container-title']) > 0:
            pubs[r]['journal'] = crossref_records[r]['container-title'][0]

        # date can show up in two different places!
        if 'year' not in pubs[r]:
            if 'published-print' in crossref_records[r]:
                year = crossref_records[r]['published-print']['date-parts'][0][0]
            elif 'journal-issue' in crossref_records[r]:
                journal_issue = crossref_records[r]['journal-issue']
                if 'published-print' in journal_issue:
                    year = journal_issue['published-print']['date-parts'][0][0]
                else:
                    year = journal_issue['published-online']['date-parts'][0][0]

            pubs[r]['year'] = int(year)

        # convert author list to pubmed format
        if 'authors' in pubs[r]:
            authors = [i.lstrip(' ') for i in pubs[r]['authors'].split(',')]
        else:
            authors = []
            for author in crossref_records[r]['author']:
                if 'given' not in author or 'family' not in author:
                    continue
                given_split = author['given'].split(' ')
                if len(given_split) > 1:
                    initials = ''.join([i[0] for i in given_split])
                else:
                    initials = given_split[0][0]
                entry = '%s %s' % (author['family'], initials)
                authors.append(entry)
        if len(authors) > etal_thresh:
            pubs[r]['authors'] = '%s et al.' % authors[0]
        else:
            pubs[r]['authors'] = ', '.join(authors)
        pubs[r]['idType'] = 'DOI'
    return(pubs)
