# Flexible Taxonomy Databases (FlexTaxD)

Flexible Taxonomy Databases - A cross platform tool for customization and merging of various taxonomic classification sources.

Supported sources in version later than v0.1.1:
* QIIME
* NCBI
* CanSNPer
* TSV

The flextaxd (flextaxd) script allows customization of databases from NCBI, QIIME or CanSNPer sources and supports export functions into NCBI formatted names and nodes.dmp files as well as a standard tab separated file (or a selected separation). The script was initially written to allow the use of GTDB with some custom modifications to allow increased resolution of selected subgroups. GTDB was created by an Australian group aimed to restructure the taxonomy relation from the NCBI taxonomy tree to strictly follow a phylogenetic structure (http://gtdb.ecogenomic.org/) this script can use the bac120_taxonomy_r89.tsv files from the GTDB downloads page as input (with the --taxonomy_type selected as QIIME). By default the script will read a Tab separated file containing parent and child (defined by column headers). The script also allows customization of the database using multiple sources and databases can be merged at a selected node(s) there is also an option to add resolution to certain subgroups (ie combine the different database types) using a tab separated file (format described below).

All data is kept in a sqlite3 database (.ftd by default) and can be dumped to NCBI formatted names and nodes.dmp files. Supported export formats are NCBI and TSV). The TSV dump format is similar to the NCBI dump except that it contains a header (parent<tab>child), has parent on the left and only uses tab to separate each column (not \<tab\>|\<tab\>).

# Installation
```
## Using pip or conda
pip3 install flextaxd
conda install flextaxd

## Using python
python setup.py install
```
# Usage
A new database(sqlite3) file will be created automatically when a taxonomy file is supplied default full path (.ftd)

Download the latest GTDB files (latest at the time of this update is "bac120_taxonomy_r89.tsv"
check https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ for latest version)
```wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac[120]_taxonomy.tsv```
#### Save the default FlexTaxD into a custom taxonomy database
```
flextaxd --taxonomy_file bac_taxonomy_r86.tsv --taxonomy_type QIIME --database .ftd
```

### Write the database into NCBI formatted nodes and names.dmp
```
flextaxd --dump
```

### Optional parameters
Use the --help option for a complete list of parameters
```
flextaxd --help
```

## Modify your database
The database update function can use either a previously built flextaxd database or directly through a TAB separated text file with headers (parent, child, (level))). Using the --parent parameter, all nodes/edges subsequent to that parent will be added (or can replace an existing node see options) with the links supplied. The parent node must exist in the database/tables and must have the same name (ex "<i>Francisella tularensis</i>"). Using the (--replace) parameter all children in the old database under the given parent will be removed, if you only want to replace for example <i>Francisella tularensis</i> be sure not to choose <i>Francisella</i> as parent.


#### Modify the database and add sub-species specifications (for <i>Francisella tularensis</i>)
```
flextaxd --mod_file custom_modification.txt --parent "Francisella tularensis" --genomeid2taxid custom_genome_annotations.txt
```
#### Example of custom_modifications.txt and custom_genome_annotations.txt
The modification file must contain three columns \<header\> and \<node\> and \<level\>) -> (Note that these tags in the file below are only there to show what is what in the file, also any number of extra tabs/spaces are there only to visualize columns! Only headers, node names and \\t \\n chars should be in the file).
```
<header>parent                          \tchild                                     \tlevel\n
<node>Francisella tularensis             \tFrancisella tularensis tularensis         \tsubspecies\n
<node>Francisella tularensis tularensis  \tB.6                                       \tsubsubspecies\n
<node>Francisella tularensis             \tFrancisella tularensis holarctica         \tsubspecies\n
```
The genome annotation file will then contain the genome_id to taxonomy(node) name as annotation. The genome id has to match the names of the genomes in the genomes_path. In particular if with use of the create_kraken_database subscript. If the genome is already annotated, the annotation will be updated.
```
GCF_00005111.1\tFrancisella tularensis tularensis
GCF_00005211.1\tFrancisella tularensis tularensis
```

## One liner version
Create, modify and dump your database
```
flextaxd --taxonomy_file bac_taxonomy_r86.tsv --taxonomy_type QIIME --mod_file custom_modification.txt --genomeid2taxid custom_genome_annotations.txt --parent "Francisella tularensis" --dump
```

#####
#    Customize the NCBI taxonomy tree
#####

The most common database to start with is the NCBI taxonomy tree, however there are many known caveats to the NCBI tree,
FlexTaxD allows modifications of the NCBI taxonomy by replacing nodes with correct structures.

Creating a custom taxonomy database using the NCBI taxonomy tree instead of FlexTaxD as base
```
Required files from NCBI (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy):
from taxdmp.zip
    names.dmp
    nodes.dmp
nucl_gb.accession2taxid.gz (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)
```
### Create a custom taxonomy database from the NCBI taxonomy
```
flextaxd
    --taxonomy_file taxonomy/nodes.dmp  ## Path to taxonomy nodes.dmp
    --taxonomy_type NCBI                ## NCBI formatted input
    --genomes_path refseq/bacteria/     ## path to ncbi-genome-download folder with bacteria (or other folder structures)
    --genomeid2taxid taxonomy/nucl_gb.accession2taxid.gz ## accession numbers to taxid annotation
    --database NCBI_taxonomy.db         ## Name of the database where all information will be kept (for future use or reuse of flextaxd)
    -o NCBI_database                    ## Output folder
    --dump                              ## Print out names.dmp and nodes.dmp from flextaxd database
```
### Modify the NCBI database using a previously created flextaxd (example database from another source (CanSNPer.db))
Replace the Francisella node in the NCBI database with the node structure from a CanSNP database containing the Francisella CanSNPer tree.
```
## Step one: create your database containing required modifications
flextaxd --taxonomy_file cansnper_tree.txt --taxonomy_type CanSNPer --genomeid2taxid cansnper_genometotaxid_annotation.txt --database canSNPer_database/CanSNPer.db

## Step two: create your custom taxonomy database from NCBI
flextaxd --database NCBI_taxonomy.db --mod_database canSNPer_database/CanSNPer.db --parent "Francisella"

## Step three: dump your database
flextaxd --dump
```

#####
#    Create a kraken database
#####

Finally there is a quick option to create a kraken2 or a krakenuniq database using your custom taxonomy database by only supplying genome names matching your annotation table given as --genomeid2taxid
Requirements: kraken2 or krakenuniq needs to be installed, For the FlexTaxD standard database, source data from genbank or refseq is required (for genomeid2taxid match)

Note: If your genome names are different, you can create a custom genome2taxid file and import into your database to match the names of your genome fasta/fna/fa files. Note that the genome files needs to be gzipped (and end with .gz) in their stored location.

### Create Kraken DB
First dump your CTDB into names.dmp and nodes.dmp (default) if that was not already done.
```
flextaxd --dump -o NCBI_database
```

Add genomes to a kraken database <my_custom_krakendb> (and create the kraken database using --create_db)
```
flextaxd-create --kraken_db my_custom_krakendb --genomes_path /path/to/genomes/folder --create_db --krakenversion kraken2
```

The script will find all fasta, fna, fa files in the given path and then add them to the krakendb, if --create_db parameter is given
the script will execute the kraken-build --build command.

### Create kraken2 NCBI database (approx 60min with 40 cores with complete bacterial genomes from genbank may 2019 (~9000))  
```
conda activate kraken2
flextaxd-create
    --database NCBI_taxonomy.db         ## flextaxd database file
    -o NCBI_database                    ## output dir (must be the same as where names and nodes.dmp were exported (--dump))
    --kraken_db NCBI_krakendb           ## Name of kraken database, absolute path not realtive to --outdir (-o)
    --genomes_path refseq/bacteria/     ## Path to ncbi-genome-download folder with bacteria (or other folder structures) using the same as when
                                        ## the database was created is recommended otw genomes may be missed
    --create_db                         ## Request script to attempt to create a kraken database (kraken2 or krakenuniq must be installed)
    --krakenversion kraken2             ## Select version of kraken which is installed
    -p 40                               ## Number of cores to use for adding genomes to the kraken database as well as the number of cores
                                            used to run kraken*-build (\*2,uniq)
    --skip "taxid"                      ## exclude genomes in taxid see details below
```

### Exclude genomes in taxid on database creation
Remove branches, the --skip parameter was implemented for benchmarking purposes as an option to remove branches by taxid, all children of the given taxid will be excluded.

flextaxd-create --skip "taxid,taxid2"
