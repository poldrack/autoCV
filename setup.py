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
VERSION = '1.3.3'


if __name__ == "__main__":

    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    setup(
        name=DISTNAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONGDESCRIPTION,
        license=LICENSE,
        version=VERSION,
        url=URL,
        download_url=DOWNLOAD_URL,
        packages=find_packages(),
        package_data={'autocv': ['templates/*', 'testdata/*']},
        scripts=[
            'bin/autoCV', 'bin/get_NSF_collaborators'],
        install_requires=[
            "pandas",
            "numpy",
            "Bio",
            "requests",
            "crossrefapi",
            "scholarly",
            "pypatent",
            "pytest"],
        classifiers=[
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.6',
            'License :: OSI Approved :: BSD License',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: MacOS'])
