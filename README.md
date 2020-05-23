# autoCV

A tool for automatic generation of a LaTeX-based curriculum vitae (CV)

### Motivation

I recently wanted to update my CV to include all of my open science activities, such as links to open access papers, code/data, and include DOIs.  Rather than doing this by hand for each publication, I decided to build an automated tool to generate a CV using PubMed and ORCID to download the publication information, and a set of text files containing other info.  It's still a work in progress but it might be helpful for you; so far I have only tested it on my own CV, and it will almost certainly need work for others (especially if you have a common name that is not uniquely identified with a simple Pubmed query). It will be most useful for more advanced researchers whose CV may be many pages long.

### Structure

The idea behind this package is to use PubMed and ORCID to obtain an up-to-date CV in a relatively automated way.
Using it requires that the user first enter some relevant information into their ORCID account:

* Education
* Employment
* Invited Positions and Distinctions
* Membership and Service

In addition, it requires generating several CSV files containing other information that is not well organized or available within ORCID:

* **[conference.csv](conference.csv)**: Conference presentations
* **[talks.csv](talks.csv)**: Colloquium and other talks
* **[funding.csv](funding.csv)**: Grants and other funding
* **[editorial.csv](editorial.csv)**: Editorial duties and reviewing
* **[additional_pubs.csv](additional_pubs.csv)**: Publications that are not found in PubMed/ORCID (including books, book chapters, and conference proceedings - note that ORCID allows addition of books but the metadata are a bit screwy, so I prefer entering them manually in this file)
* **[teaching.csv](teaching.csv)**: Courses taught

It also allows addition of links to any reference using a csv file called **[links.csv](links.csv)**.

You will need to take a look at the examples of these files in the repository to see their structure.

## Running the code

The easiest way to run the code is using Docker, which removes the need to install a full LaTeX installation.  After [installing the Docker client](https://docs.docker.com/get-docker/), you can simply use this command:

```make cv```

This will download the data and generate the CV files, and then render them to PDF.
