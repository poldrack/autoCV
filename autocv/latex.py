"""
create and render CV
"""


class LatexCV:

    def __init__(self, researcher):
        self.researcher = researcher
        assert len(researcher.publications) > 0

        self.publications = None
        self.talks = None
        self.presentations = None
        self.patents = None
        self.education = None
        self.employment = None
        self.distinctions = None
        self.service = None
        self.memberships = None
        self.teaching = None
        self.heading = None
        self.editorial = None
        self.funding = None

    def render_latex(self):
        self.render_heading()
        self.render_patents()


    def render_heading(self):
        self.heading = [
            '\\reversemarginpar \n',
            '{\LARGE %s %s. %s}\\\\[4mm] \n'
            % (  # noqa
                self.researcher.firstname.title(),
                self.researcher.middlename.title()[0],
                self.researcher.lastname.title(),
            ),
            '\\vspace{-1cm} \n\n',
            '\\begin{multicols}{2} \n',
        ]

        for a in self.researcher.address:
            self.heading.append(a + '\\\\\n')
        self.heading.append('\\columnbreak \n\n')
        self.heading.append('Phone: %s \\\\\n' % self.researcher.phone)
        self.heading.append('email: %s \\\\\n' % self.researcher.email)
        self.heading.append('url: \\href{%s}{%s} \\\\\n' % (
            self.researcher.url.replace('http://', ''), self.researcher.url))
        if self.researcher.github is not None:
            self.heading.append('url: \\href{%s}{%s} \\\\\n' % (
                self.researcher.github, self.researcher.github.replace('http://', '')))
        if self.researcher.twitter is not None:
            self.heading.append('Twitter: %s \\\\\n' % self.researcher.twitter)
        self.heading.append('ORCID: \\href{https://orcid.org/%s}{%s} \\\\\n' % (
            self.researcher.orcid, self.researcher.orcid))
        self.heading.append('\end{multicols}\n\n')  # noqa

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
