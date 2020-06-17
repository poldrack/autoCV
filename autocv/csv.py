"""
functions to read info from csv files
"""

import os
import pandas as pd
import numpy as np
import random
import string
from pprint import pprint
from collections import OrderedDict


# get code/data/osf links
def get_links_from_csv(link_file, verbose=False):
    links = {}
    if os.path.exists(link_file):
        data_df = pd.read_csv(link_file)
        for i in data_df.index:
            linktype = data_df.loc[i, 'type']
            id = data_df.loc[i, 'DOI']
            if linktype not in links:
                links[linktype] = {}
            links[linktype][id] = data_df.loc[i, 'url']
        if verbose:
            print('Links:')
            pprint(links)
    return(links)


def add_additional_pubs_from_csv(pubs, pubfile='additional_pubs.csv'):
    if os.path.exists(pubfile):
        addpubs = pd.read_csv(pubfile)
        addpubs = addpubs.fillna('')
        # resolve duplicate ISBNs (e.g. multiple chapters in a book)
        if 'ISBN' in addpubs.columns:
            isbn_loc = np.where(addpubs.columns == 'ISBN')[0][0]
            addpubs.loc[:, 'ISBN'] = [i.strip(' ') for i in addpubs.ISBN]
            for i in range(1, addpubs.shape[0]):
                if addpubs.iloc[i, isbn_loc] == '':
                    continue
                if addpubs.iloc[i, isbn_loc] in addpubs.iloc[:(i - 1), isbn_loc].tolist():
                    print('found match')
                    addpubs.iloc[i, isbn_loc] = addpubs.iloc[i, isbn_loc] + '-' + ''.join(
                        random.choice(string.ascii_lowercase) for i in range(3))
        for i in addpubs.index:
            # make a random string to stand in for pmid
            if addpubs.loc[i, 'DOI'] != '':
                id = addpubs.loc[i, 'DOI']
                idType = 'DOI'
            elif addpubs.loc[i, 'ISBN'] != '':
                id = addpubs.loc[i, 'ISBN']
                idType = 'ISBN'
            else:
                id = ''.join(random.choice(
                    string.ascii_lowercase) for i in range(8))
                idType = 'randomID'
            if id in pubs and pubs[id]['title'] == addpubs.loc[i, 'title']:
                print('found duplicate pub - skipping:', id)
                continue
            elif id in pubs:
                print('found duplicate id - modifying:', id)
                print(pubs[id])
                print('')
                id += '-' + ''.join(random.choice(
                    string.ascii_lowercase) for i in range(8))
            else:
                print('found unique id:', id)
                print('')
            pubs[id] = {idType: id}
            for c in addpubs.columns:
                if c in ['DOI', 'ISBN']:
                    continue
                entry = addpubs.loc[i, c]
                if isinstance(entry, str):
                    entry = entry.strip(' ')
                pubs[id][c] = entry
    return(pubs)


def get_teaching_from_csv(teaching_file='teaching.csv'):
    teaching_dict = OrderedDict()
    if os.path.exists(teaching_file):
        teaching_df = pd.read_csv(teaching_file)
        for i in teaching_df.index:
            coursetype = teaching_df.loc[i, 'type']
            if coursetype not in teaching_dict:
                teaching_dict[coursetype] = []
            teaching_dict[coursetype].append(teaching_df.loc[i, 'name'])

    return(teaching_dict)


def get_funding_from_csv(funding_file='funding.csv'):
    if os.path.exists(funding_file):
        grants_df = pd.read_csv(funding_file)
        return(grants_df)
    else:
        return(None)
