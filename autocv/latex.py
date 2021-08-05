"""
create and render CV
"""

import pkg_resources
import autocv.orcid as orcid
import os
import pandas as pd
from collections import OrderedDict
from autocv.utils import make_funding_line, get_links
from autocv.utils import get_pubs_by_year, get_keys_sorted_by_author, escape_characters_for_latex


class LatexCV:

    def __init__(self, researcher, font='TeX Gyre Termes', etalthresh=10, etalnum=3):
        assert len(researcher.publications) > 0
        self.researcher = researcher
        self.font = font
        self.etalthresh = etalthresh
        self.etalnum = etalnum
        self.sections_to_write = [
            'heading', 'education', 'employment', 'distinctions',
            'editorial', 'memberships', 'service',
            'funding', 'teaching', 'patents',
            'publications', 'presentations', 'talks']
        for section in self.sections_to_write:
            setattr(self, section, 'None')

        self.heading = None
        self.front_file = None
        self.back_file = None

    def load_template_files(self):
        self.front_file = pkg_resources.resource_filename(
            'autocv', 'templates/front.tex')
        with open(self.front_file) as f:
            self.front = f.readlines()
        self.back_file = pkg_resources.resource_filename(
            'autocv', 'templates/back.tex')
        with open(self.back_file) as f:
            self.back = f.readlines()

    def render_latex(self):
        for section in self.sections_to_write:
            if hasattr(self, 'render_%s' % section):
                print('rendering', section)
                render_func = getattr(self, 'render_%s' % section)
                render_func()

    def write_latex(self, outfile='autocv.tex'):
        with open(outfile, 'w') as f:
            f.write(''.join(self.front))
            for section in self.sections_to_write:
                if getattr(self, section) != 'None':
                    f.write(getattr(self, section))
            f.write(''.join(self.back))

    def render_heading(self):
        self.heading = '\\setmainfont[Ligatures=TeX]{%s} \n' \
            '\\begin{document} \n' \
            '\\reversemarginpar \n' \
            '{\\LARGE %s %s. %s}\\\\[4mm] \n' \
            '\\vspace{-1cm} \n\n' \
            '\\begin{multicols}{2} \n' % (  # noqa
                self.font,
                self.researcher.firstname.title(),
                self.researcher.middlename.title()[0],
                self.researcher.lastname.title(),
            )

        for a in self.researcher.address:
            self.heading += a + '\\\\\n'
        self.heading += '\\columnbreak \n\n'
        self.heading += 'Phone: %s \\\\\n' % self.researcher.phone
        self.heading += 'email: %s \\\\\n' % self.researcher.email
        self.heading += 'url: \\href{%s}{%s} \\\\\n' % (
            self.researcher.url.replace('http://', ''), self.researcher.url)
        if self.researcher.github is not None:
            self.heading += 'url: \\href{%s}{%s} \\\\\n' % (
                self.researcher.github, self.researcher.github.replace('http://', ''))
        if self.researcher.twitter is not None:
            self.heading += 'Twitter: %s \\\\\n' % self.researcher.twitter
        self.heading += 'ORCID: \\href{https://orcid.org/%s}{%s} \\\\\n' % (
            self.researcher.orcid, self.researcher.orcid)
        self.heading += '\\end{multicols}\n\n\\hrule\n\n'  # noqa

    def render_patents(self):
        if self.researcher.patent_data is None or len(self.researcher.patent_data) == 0:
            return
        self.patents = ('\\section*{Patents}\n\\noindent\n\n')
        for p in self.researcher.patent_data:
            authorlist = []
            for a in p['inventors']:
                initials = [i[0] for i in a[0].split(' ')]
                ln = a[1]
                authorlist.append('%s %s' % (ln, ''.join(initials)))
            authors = ', '.join(authorlist)

            self.patents += '%s (%s) \\textit{%s} US Patent \\# \\href{%s}{%s} \\vspace{2mm}\n\n' % ( # noqa
                authors,
                p['patent_date'],
                p['title'],
                p['url'],
                p['patent_num']
            )

    def render_education(self):
        education_df = orcid.get_orcid_education(self.researcher.orcid_data)
        self.education = '\\section*{Education and training}\n\\noindent\n\n'
        for i in education_df.index:
            if education_df.loc[i, 'start_date'] is None:
                continue
            self.education += '\\textit{%s-%s}: %s (%s), %s, %s\n\n' % (
                education_df.loc[i, 'start_date'],
                education_df.loc[i, 'end_date'],
                education_df.loc[i, 'degree'],
                education_df.loc[i, 'dept'],
                education_df.loc[i, 'institution'],
                education_df.loc[i, 'city'],)

    def render_employment(self):
        employment_df = orcid.get_orcid_employment(self.researcher.orcid_data)
        self.employment = '\\section*{Employment and professional affiliations}\n\\noindent\n\n'
        for i in employment_df.index:
            if employment_df.loc[i, 'dept'] is None:
                dept = ''
            else:
                dept = ' (%s)' % employment_df.loc[i, 'dept']

            self.employment += '\\textit{%s-%s}: %s%s, %s\n\n' % (
                employment_df.loc[i, 'start_date'],
                employment_df.loc[i, 'end_date'],
                employment_df.loc[i, 'role'],
                dept,
                employment_df.loc[i, 'institution'],
            )

    def render_distinctions(self):
        distinctions_df = orcid.get_orcid_distinctions(self.researcher.orcid_data)
        self.distinctions = '\\section*{Honors and Awards}\n\\noindent\n\n'
        for i in distinctions_df.index:
            self.distinctions += '\\textit{%s}: %s, %s\n\n' % (
                distinctions_df.loc[i, 'start_date'],
                distinctions_df.loc[i, 'title'],
                distinctions_df.loc[i, 'organization'],
            )

    def render_editorial(self, editorial_filename='editorial.csv'):
        editorial_file = os.path.join(self.researcher.basedir, editorial_filename)
        if os.path.exists(editorial_file):
            editorial_df = pd.read_csv(editorial_file)
        else:
            return

        editorial_df = editorial_df.fillna('')
        editorial_dict = OrderedDict()
        for i in editorial_df.index:
            role = editorial_df.loc[i, 'role']
            if role not in editorial_dict:
                editorial_dict[role] = []
            if editorial_df.loc[i, 'dates'] != '':
                date_string = ' (%s)' % editorial_df.loc[i, 'dates']
            else:
                date_string = ''
            editorial_dict[role].append(editorial_df.loc[i, 'journal'].strip(' ') + date_string)

        self.editorial = '\\section*{Editorial Duties and Reviewing} \n\\noindent \n\n'
        for i in editorial_dict:
            self.editorial += '\\textit{%s}: %s \n\n' % (
                i.strip(' '), ', '.join(editorial_dict[i]))

    def render_service(self):
        service_df = orcid.get_orcid_service(self.researcher.orcid_data)

        self.service = '\\section*{Service}\n\\noindent\n\n'
        for i in service_df.index:
            self.service += '%s, %s, %s-%s \n\n' % (
                service_df.loc[i, 'role'],
                service_df.loc[i, 'organization'],
                service_df.loc[i, 'start_date'],
                service_df.loc[i, 'end_date'],
            )

    def render_memberships(self):
        memberships_df = orcid.get_orcid_memberships(self.researcher.orcid_data)
        self.memberships = '\\section*{Professional societies}\n\\noindent\n\n'
        self.memberships += ', '.join(memberships_df.organization) + '\n\n'

    def render_teaching(self, teaching_filename='teaching.csv'):
        teaching_file = os.path.join(self.researcher.basedir, teaching_filename)
        if os.path.exists(teaching_file):
            teaching_df = pd.read_csv(teaching_file)
        else:
            return

        teaching_dict = OrderedDict()
        for i in teaching_df.index:
            coursetype = teaching_df.loc[i, 'type']
            if coursetype not in teaching_dict:
                teaching_dict[coursetype] = []
            teaching_dict[coursetype].append(teaching_df.loc[i, 'name'])

        self.teaching = '\\section*{Teaching}\n\\noindent\n\n'
        for i in teaching_dict:
            self.teaching += '\\textit{%s}: %s \\vspace{2mm}\n\n' % (i, ', '.join(teaching_dict[i]))

    def render_funding(self, funding_filename='funding.csv', abbreviate=True):
        funding_file = os.path.join(self.researcher.basedir, funding_filename)
        if os.path.exists(funding_file):
            funding_df = pd.read_csv(funding_file).fillna('')
        else:
            return

        self.funding = '\\section*{Research funding}\n\\noindent\n\n'

        self.funding += '\\subsection*{Active:}\n\n'

        active_pi_funding = funding_df.query('active == True and role == "Principal Investigator"')
        for i in active_pi_funding.index:
            line = make_funding_line(active_pi_funding, i, abbreviate)
            self.funding += '%s \\vspace{2mm}\n\n' % line

        active_coi_funding = funding_df.query('active == True and role != "Principal Investigator"')
        for i in active_coi_funding.index:
            line = make_funding_line(active_coi_funding, i, abbreviate)
            self.funding += '%s \\vspace{2mm}\n\n' % line

        self.funding += '\\subsection*{Completed:}\n\n'
        completed_funding = funding_df.query('active == False')
        for i in completed_funding.index:
            line = make_funding_line(completed_funding, i, abbreviate)
            self.funding += '%s \\vspace{2mm}\n\n' % line

    def render_publications(self):
        linkfile = os.path.join(self.researcher.basedir, 'links.csv')
        links = get_links(linkfile)

        years = list({self.researcher.publications[i].year for i in self.researcher.publications})
        years.sort(reverse=True)
        if hasattr(self.researcher, 'gscholar'):
            hindex = self.researcher.gscholar.hindex
        elif hasattr(self.researcher, 'gscholar_data'):
            # loaded from json
            try:
                hindex = self.researcher.gscholar_data['hindex']
            except TypeError:
                hindex = self.researcher.gscholar_data.hindex

        self.publications = '\\section*{Publications (Google Scholar H-index = %d)}' % hindex

        for year in years:
            self.publications += '\\subsection*{%s}' % year

            year_pubs = get_pubs_by_year(self.researcher.publications, year)

            keys_sorted_by_author = get_keys_sorted_by_author(year_pubs)

            # get the reference line for each
            for pub in keys_sorted_by_author:
                self.researcher.publications[pub] = escape_characters_for_latex(self.researcher.publications[pub])

                line = self.researcher.publications[pub].format_reference_latex(self.etalthresh, self.etalnum)

                if hasattr(self.researcher.publications[pub], 'PMC') and self.researcher.publications[pub].PMC is not None:
                    line += ' \\href{https://www.ncbi.nlm.nih.gov/pmc/articles/%s}{OA}' % self.researcher.publications[pub].PMC

                for linktype in links:
                    if pub in links[linktype]:
                        line += ' \\href{%s}{%s}' % (
                            links[linktype][pub], linktype)

                # TBD: need to filter out non-DOIs
                if self.researcher.publications[pub].type not in ['book', 'monograph'] and hasattr(self.researcher.publications[pub], 'DOI'):
                    line += ' \\href{http://dx.doi.org/%s}{DOI}' % self.researcher.publications[pub].DOI

                line += ' \\vspace{2mm}\n\n'
                self.publications += line

    def render_presentations(self, presentations_filename='conference.csv'):
        presentations_file = os.path.join(self.researcher.basedir, presentations_filename)
        if os.path.exists(presentations_file):
            presentations = pd.read_csv(presentations_file) #, index_col=0)
        else:
            return

        presentations = presentations.sort_values('year', ascending=False)

        self.presentations = '\\section*{Conference Presentations}\n\\noindent\n\n'

        for i in presentations.index:
            entry = presentations.loc[i, :]
            title = entry.title.strip('.')
            location = entry.location.strip(' ').strip('.')
            self.presentations += '%s (%s). \\emph{%s}. %s. \\vspace{2mm}\n\n' % (
                entry.authors, entry.year, title, location)

    def render_talks(self, talks_filename='talks.csv', verbose=True):
        talks_file = os.path.join(self.researcher.basedir, talks_filename)
        if os.path.exists(talks_file):
            talks = pd.read_csv(talks_file) #, index_col=0)
        else:
            return
        years = list(talks.year.unique())
        years.sort()
        if verbose:
            years = years[::-1]
        print('years:', years)
        lines = '\\section*{Invited addresses and colloquia (* - talks given virtually)}\n\\noindent\n\n'
        for y in years:
            if verbose:
                print('year', y)
            talks_year = talks.query('year == %s' % y)
            lines += '%s: %s \\vspace{2mm}\n\n' % (y, ','.join(list(talks_year.place)))
        self.talks = lines
