#!/usr/bin/env python

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

def get_version():
    with open(path.join(here, "centr_utils/version.py"), encoding = 'utf-8') as hin:
        for line in hin:
            if line.startswith("__version__"):
                version = line.partition('=')[2]
                return version.strip().strip('\'"')
    raise ValueError('Could not find version.')

setup(
    name = 'centr_utils',
    version = get_version(),
    description='Python tools for analyzing centromere sequences',
    long_description=long_description, 
    long_description_content_type='text/markdown',
    url = 'https://github.com/friend1ws/centr_utils',
    author = 'Yuichi Shiraishi',
    author_email = 'friend1ws@gamil.com',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['tests']),
    install_requires = [],

    # packages = find_packages(exclude = ['tests', 'LINE_db']),
    # package_data={'nanomonsv': ['data/*']},

    # install_requires = ["numpy", "parasail", "pysam"],
    entry_points = {'console_scripts': ['centr_utils = centr_utils:main']}

)


