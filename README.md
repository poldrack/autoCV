# autoCV
![Python package](https://github.com/poldrack/autoCV/workflows/Python%20package/badge.svg)
[![codecov](https://codecov.io/gh/poldrack/autoCV/branch/master/graph/badge.svg)](https://codecov.io/gh/poldrack/autoCV)
[![PyPI version](https://badge.fury.io/py/autocv.svg)](https://badge.fury.io/py/autocv)

A tool for automatic generation of a LaTeX-based curriculum vitae (CV)

### Motivation

I recently wanted to update my CV to include all of my open science activities, such as links to open access papers, code/data, and include DOIs.  Rather than doing this by hand for each publication, I decided to build an automated tool to generate a CV using PubMed and ORCID to download the publication information, and a set of text files containing other info.  It's still a work in progress but it might be helpful for you; so far I have only tested it on my own CV, and it will almost certainly need work for others (especially if you have a common name that is not uniquely identified with a simple Pubmed query). It will be most useful for more advanced researchers whose CV may be many pages long.  

The project takes advantage of the [very nice LaTeX CV template](http://nitens.org/taraborelli/cvtex) from Dario Taraborelli.

### Structure

The idea behind this package is to use PubMed and ORCID to obtain an up-to-date CV in a relatively automated way.
Using it requires that the user first enter some relevant information into their ORCID account:

* Education
* Employment
* Invited Positions and Distinctions
* Membership and Service

In addition, it requires generating several CSV files containing other information that is not well organized or available within ORCID:

* **[conference.csv](autocv/testdata/conference.csv)**: Conference presentations
* **[talks.csv](autocv/testdata/talks.csv)**: Colloquium and other talks
* **[funding.csv](autocv/testdata/funding.csv)**: Grants and other funding
* **[editorial.csv](autocv/testdata/editorial.csv)**: Editorial duties and reviewing
* **[additional_pubs.csv](autocv/testdata/additional_pubs.csv)**: Publications that are not found in PubMed/ORCID (including books, book chapters, and conference proceedings - note that ORCID allows addition of books but the metadata are a bit screwy, so I prefer entering them manually in this file)
* **[teaching.csv](autocv/testdata/teaching.csv)**: Courses taught

It also allows addition of links to any reference using a csv file called **[links.csv](autocv/testdata/links.csv)**.

Finally, you will need to generate a json file called params.json that contains some metadata about you - see example [here](autocv/testdata/params.json).

You will need to take a look at the examples of these files in the repository to see their structure.

You can see an example of the output [here](autocv/testdata/autocv_template.pdf).

## Running the code

First you should install the package:

```pip install -U autocv```

Then you can run the full process by simply typing ```autoCV``` from the directory that contains all of the necessary files.  Type ```autoCV -h``` for additional options.