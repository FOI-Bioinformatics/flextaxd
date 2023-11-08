# Flexible Taxonomy Databases (FlexTaxD)

FlexTaxD is a versatile tool for the customization and integration of taxonomic classifications from diverse sources. It facilitates the creation, modification, and export of taxonomy databases for bioinformatics applications.

## Key Features
- **Supported Taxonomy Formats**: NCBI, TSV, MP-style (GTDB, QIIME, SILVA, CanSNPer).
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
