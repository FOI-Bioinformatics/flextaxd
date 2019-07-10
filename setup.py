from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from custom_taxonomy.custom_taxonomy_databases import __version__

import sys
if sys.version_info.major < 3 and sys.version_info.minor < 5:
    current_version = ".".join(map(str,[sys.version_info.major,sys.version_info.minor,sys.version_info.micro]))
    exit("This script only supports python versions 3.5 and above, please upgrade python! Current version: {python} ".format(python=current_version ))

localpath = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(localpath, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='custom_taxonomy_databases',
    ##Global version, does not nessesarily follow script versions
    version=__version__,

    description='Script that allows the creation of custom kraken databases from various sources (NCBI, QIIME, CanSNPer)',
    long_description=long_description,
    long_description_content_type="text/markdown",

    # The project's main homepage.
    url='https://github.com/davve2/custom-taxonomy-databases',

    # Author details
    author='David Sundell',
    author_email='david.sundell@foi.se',

    # Choose your license
    license=' GNU GENERAL PUBLIC LICENSE version 3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: Beta',

        'Intended Audience :: Developers bioinformatics',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved ::  GNU GENERAL PUBLIC LICENSE version 3 License',

        'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7'
    ],

    keywords='taxonomy NCBI CanSNPer GTDB customization',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    #packages=find_packages(exclude=['contrib', 'docs', 'test*']),
	packages=find_packages(exclude=['contrib', 'docs', 'test*']),
	#py_modules=["custom_taxonomy.create_taxonomy_databases"]
	#packages=['source','modules'],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'custom_taxonomy_databases=custom_taxonomy.custom_taxonomy_databases:main',
			'create_databases=custom_taxonomy.create_databases:main',
        ],
    },
)
