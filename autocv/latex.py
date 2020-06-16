"""
create and render CV
"""

import pkg_resources
import autocv.orcid as orcid


class LatexCV:

    def __init__(self, researcher, font='TeX Gyre Termes'):
        assert len(researcher.publications) > 0
        self.researcher = researcher
        self.font = font
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
            '{\LARGE %s %s. %s}\\\\[4mm] \n' \
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
        self.heading += '\end{multicols}\n\n\hrule\n\n'  # noqa

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

            self.patents += '%s (%s) \\textit{%s} US Patent \# \\href{%s}{%s} \\vspace{2mm}\n\n' % ( # noqa
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
