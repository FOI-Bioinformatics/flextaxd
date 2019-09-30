from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from flextaxd.custom_taxonomy_databases import __version__

import sys
if sys.version_info.major < 3 and sys.version_info.minor < 5:
    current_version = ".".join(map(str,[sys.version_info.major,sys.version_info.minor,sys.version_info.micro]))
    exit("This script only supports python versions 3.5 and above, please upgrade python! Current version: {python} ".format(python=current_version ))

localpath = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(localpath, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='flextaxd',
    ##Global version, does not nessesarily follow script versions
    version=__version__,

    description='Script that allows the creation of custom kraken databases from various sources (NCBI, QIIME, CanSNPer)',
    long_description=long_description,
    long_description_content_type="text/markdown",

    # The project's main homepage.
    url='https://github.com/FOI-Bioinformatics/flextaxd',

    # Author details
    author='David Sundell',
    author_email='david.sundell@foi.se',

    # Choose your license
    license=' GNU GENERAL PUBLIC LICENSE version 3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Programming Language :: Python :: 3',

        #'Intended Audience :: Developers bioinformatics',
        #'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.6',
    keywords='taxonomy NCBI CanSNPer customization GTDB',
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().,
	packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    entry_points={
        'console_scripts': [
            'flextaxd=flextaxd.custom_taxonomy_databases:main',
			'flextaxd-create=flextaxd.create_databases:main',
        ],
    },
)
