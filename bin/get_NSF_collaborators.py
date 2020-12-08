#!/usr/bin/env python3
"""
get NSF collaborators list using autoCV
"""

import os
import pandas as pd
from collections import defaultdict
import argparse

from autocv.researcher import Researcher
from datetime import datetime
from dateutil.relativedelta import relativedelta


def convert_pubmed_date_to_datetime(pubmed_record):
    try:
        pubmed_date = pubmed_record[
            'MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
    except KeyError:  # some records don't have pub date - skip them
        return(None)

    year = int(pubmed_date['Year']) if 'Year' in pubmed_date else None
    if 'Month' not in pubmed_date:
        month = 12
    elif pubmed_date['Month'] in monthnums:
        month = monthnums[pubmed_date['Month']] if 'Month' in pubmed_date else None
    else:
        month = int(pubmed_date['Month'])

    if year is None and month is None:
        if 'MedlineDate' in pubmed_date:
            year = int(pubmed_date['MedlineDate'].split(' ')[0])
        else:
            year = None
    retval = datetime(year, month, 1) if year is not None and month is not None else None
    return(retval)


def get_args():
    parser = argparse.ArgumentParser(description='NSF collaborator list generator')
    parser.add_argument('--paramfile', default='params.json',
                        help='json file containing parameters')
    parser.add_argument('--basedir', default='.',
                        help='directory for output')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    basedir = args.basedir

    assert os.path.exists(args.paramfile), "You must first set up your params.json file, see README.md for more"

    r = Researcher(args.paramfile)
    r.get_orcid_data()
    r.get_orcid_dois()
    r.get_pubmed_data()
    # r.get_google_scholar_record()
    # r.get_patents()
    # r.make_publication_records()

    # NSF requires all relationships within last 48 months
    cutoff = datetime.now() - relativedelta(years=4)
    monthnums = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
                 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}

    # first get full list of authors from pubmed articles
    skipped_bad_date = []
    coauthors = []
    for p in r.pubmed_data['PubmedArticle']:

        pubdate = convert_pubmed_date_to_datetime(p)
        if pubdate is None:
            skipped_bad_date.append(p)
            continue

        if pubdate < cutoff:  # exclude papers more than 4 years old
            continue

        pmid = int(p['MedlineCitation']['PMID'])

        if 'AuthorList' in p['MedlineCitation']['Article']:
            for au in p['MedlineCitation']['Article']['AuthorList']:
                coauthors.append((au, pubdate, pmid))

    print(f'skipped {len(skipped_bad_date)} pubs due to bad records')

    # create a dictionary keyed by author names
    authordict = {}
    skipped_bad_name = []
    first_initials = []

    for au, date, _ in coauthors:
        try:
            fullname = au['LastName'] + '_' + au['ForeName']
            initials_only = True
            for slot in au['ForeName'].split(' '):
                if len(slot) > 1:
                    initials_only = False
            if initials_only:
                first_initials.append(fullname)
        except KeyError:
            skipped_bad_name.append(au)
            continue

        dt = pd.to_datetime(date, format='%Y-%m')

        if fullname in authordict:
            authordict[fullname]['dates'].append(dt)
        else:
            # check for email address
            if len(au['AffiliationInfo']):
                affil = au['AffiliationInfo'][0]['Affiliation']
                if affil.split(' ')[-1].find('@') > -1:
                    # email=affil.split(' ')[-1]
                    affil = ' '.join(affil.split(' ')[:-1])
            else:
                affil = []
            email = ''
            authordict[fullname] = {'Affiliation': affil, 'email': email,
                                    'dates': [dt]}

    # find potential duplicates
    author_name_dict = defaultdict(lambda: [])
    for k in authordict.keys():
        lastname, firstplus = k.split('_')
        first_split = firstplus.split(' ')
        author_name_dict[lastname].append(first_split)

    # Try to resolve dupes
    for k, forenames in author_name_dict.items():
        if len(forenames) < 2:
            continue

        firstinitials = [i[0][0] for i in forenames]
        # continued if no duplicate first intials
        if len(set(firstinitials)) == len(firstinitials):
            continue
        print(k, forenames)

        forename_lengths = [len(i) for i in forenames]  # length of first + initials
        first_name_lengths = [len(i[0]) for i in forenames]  # length of first only
        middle_initial_lengths = [len(i[1]) for i in forenames if len(i) > 1]  # length of middle only

        replace = False

        # find dupes where first name is same but one has initial and the other doesn't
        if set(forename_lengths) == {1, 2} and len(forenames) == 2 and min(middle_initial_lengths) == 1:
            all_initial_idx = forename_lengths.index(min(forename_lengths))
            all_initial_key = f'{k}_{" ".join(forenames[all_initial_idx])}'
            full_name_idx = forename_lengths.index(max(forename_lengths))
            full_name_key = f'{k}_{" ".join(forenames[full_name_idx])}'
            replace = True

        # find dupes for cases where there is just one name and one initial
        # we can assume they will match due to test above
        elif set(forename_lengths) == {2} and len(forenames) == 2 and min(first_name_lengths) == 1:
            all_initial_idx = first_name_lengths.index(1)
            all_initial_key = f'{k}_{" ".join(forenames[all_initial_idx])}'
            full_name_idx = first_name_lengths.index(max(first_name_lengths))
            full_name_key = f'{k}_{" ".join(forenames[full_name_idx])}'
            replace = True

        # find dupes with two two-length entries
        # we currently punt if there are more than two names
        elif set(forename_lengths) == {2} and len(forenames) == 2:

            # make sure that both intials match
            initials = ['', '']
            for i in range(2):
                for j in range(2):
                    initials[i] += forenames[i][j][0]
            if initials[0] != initials[1]:
                # initials don't match
                continue

            # find if only one of them is all initials
            all_initials = [None, None]
            for i in range(2):
                all_initials[i] = sum([len(j) for j in forenames[i]]) == 2

            if sum(all_initials) == 1:
                all_initial_idx = all_initials.index(True)
                all_initial_key = f'{k}_{" ".join(forenames[all_initial_idx])}'
                full_name_idx = all_initials.index(False)
                full_name_key = f'{k}_{" ".join(forenames[full_name_idx])}'
                replace = True

        if replace:
            print('all initials, consolidating and removing:', forenames[all_initial_idx])
            authordict[full_name_key]['dates'] += authordict[all_initial_key]['dates']
            if authordict[full_name_key]['Affiliation'] == '':
                authordict[full_name_key]['Affiliation'] = authordict[all_initial_key]['Affiliation']
            del authordict[all_initial_key]
            # move the items from this entry to the other one,
            # and then remove this one

            print()

    # find latest date for each coauthor
    for k in authordict.keys():
        if len(authordict[k]['dates']) > 1:
            authordict[k]['dates'].sort()
        authordict[k]['latest'] = authordict[k]['dates'][-1]

    # get sorted list of names
    coauthor_names = list(authordict.keys())
    coauthor_names.sort()

    with open('collabs.txt', 'w') as f:
        for a in coauthor_names:
            f.write("A:\t%s\t\t%s\t%d/1/%d\n" % (
                a.replace('_', ', '), authordict[a]['Affiliation'],
                authordict[a]['latest'].month, authordict[a]['latest'].year))
