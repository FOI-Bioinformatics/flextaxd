# Flexible Taxonomy Databases (FlexTaxD)

FlexTaxD is a versatile tool for the customization and integration of taxonomic classifications from diverse sources. It facilitates the creation, modification, and export of taxonomy databases for bioinformatics applications.

## Key Features
- **Supported Taxonomy Formats**: NCBI, TSV, MP-style (GTDB, SILVA, CanSNPer). The lineage format MP names here (`d__Bacteria;p__...;g__...;s__...`, one genome per row) was called QIIME until v0.6.1: the file and the reader are unchanged, only the name of the option (`-tt MP`, and `-tt GTDB` still maps to it).
- **Database Build Programs**: Compatible with kraken2, ganon, krakenuniq, and centrifuge.
- **Database Customization**: Modify, annotate, and clean up taxonomy databases to fit specific research needs.
- **Output**: Export databases to NCBI formatted files or tab-separated values, with customizable options.
- **Data Management**: Utilizes a SQLite database to manage data efficiently.

## Quick Start
For a complete walkthrough, refer to the [FlexTaxD Wiki](https://github.com/FOI-Bioinformatics/flextaxd/wiki).

### Installation
```bash
# With conda using mamba
conda install mamba -n base -c conda-forge
mamba create -c conda-forge -c bioconda -n flextaxd flextaxd

# Manual Python installation
git clone https://github.com/FOI-Bioinformatics/flextaxd
cd flextaxd
pip install .
```

### Usage Examples
```bash
# Create a custom taxonomy database
flextaxd --taxonomy_file taxonomy.tsv --database .ftd

# Export database to NCBI format
flextaxd --dump

# Print database statistics
flextaxd --stats

# Get help
flextaxd --help
```

### Taxonomy dumps and the format a program reads
```bash
# Export the database in the layout kraken2/krakenuniq read
flextaxd --dump --dbprogram kraken2

# Minimal, tab separated export (no program will read this one)
flextaxd --dump_mini
```

kraken2 and krakenuniq take the name of a taxid from the trailing columns of `names.dmp`, and the rank
from the trailing columns of `nodes.dmp`. `--dump --dbprogram kraken2` writes those columns; `--dump_mini`
writes a short tab separated file without them, which is meant for reading by a person or another tool.
The two are never combined: `--dump_mini --dbprogram kraken2` (and an NCBI type dump for a program) stops
with an error instead of writing a dump that leaves the built database listing every taxid without a name.

## The two levels: taxonomy database and classifier database

`flextaxd` works on one **taxonomy database**, the sqlite file named by `-db`. It is written once and then
read back: modified with `-mf`/`-gt`, merged with `-md`, cleaned with `--clean_database`, inspected with
`--stats`/`--report` and dumped with `--dump`. A **classifier database** (kraken2, krakenuniq, ganon,
centrifuge) is the output of the second command, `flextaxd-create`, in the directory named by `--db_name`.
It is output only: a built kraken2 database cannot be imported into flextaxd. What can be imported are the
inputs of such a build — the taxdump files, a genome folder, an `accession2taxid` annotation file, and
packed sequence files given with `-nt`/`--nt_source` and `--multifile_prefix`.

`-db` is the taxonomy database in both commands, `--db_name` is the classifier database the build writes;
they are never the same path.

### From a taxonomy to a classifier database, end to end

```bash
# 1 import a taxonomy once, into the taxonomy database
flextaxd -db life.ftd -tf bac120_taxonomy.tsv -tt MP

# 2 merge a branch of another taxonomy database into it, then annotate genomes and clean the tree
flextaxd -db life.ftd -md ncbi.ftd --parent root -gt genome2node.tsv
flextaxd -db life.ftd -mf changes.tsv -gt genome2node.tsv --clean_database

# 3 see what the database and the sequence files agree on, and dump the taxonomy
flextaxd -db life.ftd --genomes_path genomes --stats --report
flextaxd -db life.ftd --dump --dbprogram kraken2 -o dumps

# 4 build the classifier database from it
flextaxd-create -db life.ftd --db_name kraken_run --dbprogram kraken2 --genomes_path genomes -nt nt.fna.gz --create_db -o build
```

Step 2 is the merge people look for when they ask to merge two kraken2 databases: merge the taxonomy
databases, then build the classifier database again. The step that decides whether a clean drops genomes
is the annotation, `--report` names which annotated genomes have no sequence file and which library
sequences have no taxid; for an NCBI source the second number is usually an `accession2taxid` file that
covers only `nucl_gb` and not `nucl_wgs`.


## Requirements
- Python >=3.6
- Additional dependencies vary based on executed functions (e.g., `ncbi-genome-download` for genome downloads).

### Visualization Dependencies
- **biopython**: Required for Newick visualizations.
- **matplotlib**: For tree visualizations.
- **inquirer**: For interactive prompts when multiple parents are present.

### Database Creation Dependencies
- kraken2, krakenuniq, ganon, centrifuge: Required if `create_database` is used.

## Contributing
Your contributions are welcome! Please refer to the [Contribution Guide](https://github.com/FOI-Bioinformatics/flextaxd/CONTRIBUTING.md) for details on how to submit pull requests, report issues, or request features.

## License
FlexTaxD is open-sourced under the [MIT license](https://github.com/FOI-Bioinformatics/flextaxd/LICENSE).



## Citation
FlexTaxD is published in [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab621/6361544).
```
Sundell, D. et al. (2021) ‘FlexTaxD: flexible modification of taxonomy databases for improved sequence classification’, Bioinformatics. Edited by J. Kelso. Bioinformatics. doi: 10.1093/bioinformatics/btab621.
```
