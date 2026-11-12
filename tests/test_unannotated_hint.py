"""issue #71: what the tool says when a library loses many of its sequences to a missing taxid

The question in the thread was why a clean dropped so many annotated genomes and why the kraken2 database
grew so much after the annotation was widened. The answer is the annotation file: NCBI describes the
assemblies of its standard set in nucl_gb.accession2taxid and the contig level assemblies and metagenome
assembled genomes in nucl_wgs.accession2taxid, so a database imported from one file places less of its
library than it could. The tool now says that when the numbers show the case, and says nothing otherwise."""

from flextaxd.modules.functions import annotation_hint


def test_the_hint_names_the_other_annotation_file_and_the_numbers():
    """the measured numbers of the issue thread, 25514 sequences of 41170 without a taxid"""
    hint = annotation_hint(41170, 25514)
    assert "25514 of 41170" in hint and "62%" in hint, hint
    assert "nucl_gb.accession2taxid" in hint and "nucl_wgs.accession2taxid" in hint, hint


def test_a_library_that_is_mostly_annotated_gets_no_advice():
    """the advice is only for the case where a large share of the library is missing, a stray sequence or
    ten per cent of a large library is not that case"""
    assert annotation_hint(41170, 20) == ""
    assert annotation_hint(41170, 4000) == ""
    assert annotation_hint(20, 20) == ""
    assert annotation_hint(0, 0) == ""


def test_the_share_and_the_count_both_have_to_be_large():
    """a quarter of the library is the share where the annotation file is the likely answer, but only from
    a count that is not a coincidence of a small file"""
    assert annotation_hint(100, 25), "a quarter of a hundred sequences is the case"
    assert annotation_hint(1000, 249) == "", "24.9 per cent is not"
    assert annotation_hint(40, 24) == "", "too few sequences to name a cause"
    assert annotation_hint(40, 25), "25 of 40 is a quarter of the library either way"
