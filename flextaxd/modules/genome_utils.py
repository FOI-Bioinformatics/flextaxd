"""Shared helpers for extracting genome IDs from filenames."""

import re

FASTA_EXT = (".fna", ".fa", ".fasta")

# Pattern 1: GCF_/GCA_ + 9 digits + optional .version
_GCX_PATTERN = r'GC[AF]_\d{9}(?:\.\d+)?'
# Pattern 2: 2 uppercase letters + 6 digits + optional .version
_ACC_PATTERN = r'[A-Z]{2}\d{6}(?:\.\d+)?'

_GENOME_ID_RE = re.compile(
    r'(?:^|_)(' + _GCX_PATTERN + r'|' + _ACC_PATTERN + r')(?:_|$)'
)


def strip_fasta_ext(filename):
    """Remove .gz suffix, then FASTA extension from filename."""
    if filename.endswith(".gz"):
        filename = filename[:-3]
    for ext in FASTA_EXT:
        if filename.endswith(ext):
            return filename[:-len(ext)]
    return filename


def extract_genome_id(filename):
    """Extract genome ID from filename (after stripping extensions).

    Returns the genome ID string, or None if no pattern matches.
    """
    name = strip_fasta_ext(filename)
    m = _GENOME_ID_RE.search(name)
    if m:
        return m.group(1)
    return None
