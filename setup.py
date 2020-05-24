#! /usr/bin/env python
#
# Copyright (C) 2013-2015 Russell Poldrack <poldrack@stanford.edu>
# some portions borrowed from https://github.com/mwaskom/lyman/blob/master/setup.py
import os
from setuptools import setup, find_packages

descr = """autoCV: An automated CV generator"""

DISTNAME = "autocv"
DESCRIPTION = descr
LONGDESCRIPTION = "This package automatically generates an academic curriculum vita (CV) using a set of open APIs along with a set of local files."
MAINTAINER = 'Russ Poldrack'
MAINTAINER_EMAIL = 'poldrack@stanford.edu'
LICENSE = 'MIT'
URL = 'https://github.com/poldrack/autoCV'
DOWNLOAD_URL = 'https://github.com/poldrack/autoCV'
VERSION = '1.0a2'


def check_dependencies():

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    needed_deps = [
        "pandas",
        "numpy",
        "Bio",
        "requests",
        "crossref",
        "scholarly",
        "pypatent",
        "pytest"]
    missing_deps = []
    for dep in needed_deps:
        try:
            __import__(dep)
        except ImportError:
            missing_deps.append(dep)

    if missing_deps:
        missing = (", ".join(missing_deps))
        raise ImportError("Missing dependencies: %s" % missing)


if __name__ == "__main__":

    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    import sys
    if not (len(sys.argv) >= 2 and (
        '--help' in sys.argv[1:] or sys.argv[1] in (
            '--help-commands',
            '--version',
            'egg_info',
            'clean'))):
        check_dependencies()

    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONGDESCRIPTION,
          license=LICENSE,
          version=VERSION,
          url=URL,
          download_url=DOWNLOAD_URL,
          packages=find_packages(),
          package_data={'autocv': ['templates/*']},
          scripts=[
              'bin/make_cv.py',
              'bin/autoCV'],
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 3.6',
              'License :: OSI Approved :: BSD License',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'])
