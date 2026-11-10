import os
import sqlite3

import pytest

from flextaxd.modules.ReadTaxonomy import ReadTaxonomy
from flextaxd.modules.ProcessDirectory import ProcessDirectory
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions, flextaxd_version


def make_annotated_database(tmp_path, genomes):
    """Create a database with a root, one child per genome and an annotation for each genome"""
    db = str(tmp_path / "report.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    taxids = {}
    for number, genome in enumerate(genomes, start=2):
        conn.query("INSERT INTO nodes (id,name) VALUES ({},?)".format(number), insert_val=(genome,))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,{},1)".format(number))
        conn.query("INSERT INTO genomes (id,genome) VALUES (?,?)", insert_val=(number, genome))
        taxids[genome] = number
    conn.commit()
    return db, taxids


def make_old_database(tmp_path):
    """A database written before flextaxd stored build information: no meta and report tables"""
    db = str(tmp_path / "old.fdb")
    conn = sqlite3.connect(db)
    conn.executescript(
        """CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
        CREATE TABLE tree (parent integer NOT NULL, child integer NOT NULL, rank_i integer, unique (parent,child,rank_i));
        CREATE TABLE rank (rank_i integer PRIMARY KEY, rank VARCHAR(15));
        CREATE TABLE genomes (id integer NOT NULL, genome VARCHAR(255) NOT NULL UNIQUE, reference VARCHAR(255));
        INSERT INTO nodes (id,name) VALUES (1,'root');
        INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1);"""
    )
    conn.commit()
    conn.close()
    return db


def test_new_database_records_the_flextaxd_version(tmp_path):
    db = str(tmp_path / "meta.fdb")
    conn = ModifyFunctions(database=db)
    assert conn.get_meta("built_with") == flextaxd_version()
    assert conn.get_meta("built")
    assert conn.get_meta("modified_with") is None
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.commit()
    assert ModifyFunctions(database=db).get_meta("modified_with") == flextaxd_version()


def test_database_without_build_information_is_upgraded_on_open(tmp_path):
    db = make_old_database(tmp_path)
    conn = ModifyFunctions(database=db)
    assert conn.get_meta("built_with") is None
    assert conn.report_kinds() == {}
    report = ModifyFunctions(database=db).report_lines()
    assert "earlier flextaxd version" in "\n".join(report)


def test_database_that_cannot_be_written_is_still_readable(tmp_path):
    """A read only database from an earlier flextaxd cannot record the build information, but
    reporting on it must still work"""
    db = make_old_database(tmp_path)
    os.chmod(db, 0o444)
    try:
        conn = ModifyFunctions(database=db)
        assert conn.get_meta("built_with") is None
        assert conn.report_lines()[0].startswith("Database built with flextaxd: unknown")
        assert conn.num_rows("nodes") == 1
    finally:
        os.chmod(db, 0o644)


def test_disk_report_lists_missing_genomes_and_unannotated_files(tmp_path):
    db, _ = make_annotated_database(tmp_path, ["GCF_000000001.1", "GCF_000000002.1"])
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_000000001.1_genomic.fna").write_text(">a\nACGT\n")
    (genomes / "GCF_000000003.1_genomic.fna").write_text(">b\nACGT\n")
    report = ProcessDirectory(db).report_disk(str(genomes))
    assert report["missing"] == ["GCF_000000002.1"]
    assert report["unannotated"] == ["GCF_000000003.1_genomic.fna"]
    assert report["annotated"] == 2 and report["on_disk"] == 1
    stored = ModifyFunctions(database=db)
    assert stored.report_kinds() == {"missing_on_disk": 1, "unannotated_file": 1}
    assert [entry for source, entry in stored.get_report("missing_on_disk")] == ["total\t1"]
    assert stored.get_meta("genomes_missing_on_disk") == "1"
    assert stored.get_meta("files_unannotated") == "1"
    lines = "\n".join(stored.report_lines())
    assert "Genome annotations without sequence file: 1" in lines
    assert "Sequence files without genome annotation: 1" in lines
    ## the database keeps the numbers, not the individual genomes and files
    assert "GCF_000000002.1" not in lines and "GCF_000000003.1_genomic.fna" not in lines


def test_sequence_level_annotations_are_not_missing_genomes(tmp_path):
    """issue #131: the packed --nt_source route annotates sequences, and by the decision that keeps those
    files in one piece they have no sequence file of their own in a genome folder. Counting them as
    annotations without a file made --report name sequences the library had in fact got from the packed
    file, made --download offer accessions for download, and made --purge_database delete the annotation
    of sequences that were used."""
    db, _ = make_annotated_database(tmp_path, ["GCF_000000001.1", "GCF_000000002.1",
                                               "NC_000913.3", "CP000250.1"])
    conn = ModifyFunctions(database=db)
    for accession in ("NC_000913.3", "CP000250.1"):
        conn.query("UPDATE genomes SET reference='nt' WHERE genome = ?", insert_val=(accession,))
    conn.commit()

    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_000000001.1_genomic.fna").write_text(">a\nACGT\n")

    _, download = ProcessDirectory(db).process_folder(str(genomes))
    assert [entry["genome_id"] for entry in download] == ["GCF_000000002.1"], \
        "the accessions of a packed file were offered for download"

    report = ProcessDirectory(db).report_disk(str(genomes))
    assert report["missing"] == ["GCF_000000002.1"]
    assert report["sequence_annotations"] == 2, "the sequence level annotations are not counted anywhere"

    stored = ModifyFunctions(database=db)
    assert stored.get_meta("genomes_missing_on_disk") == "1"
    assert stored.get_meta("sequence_annotations_on_disk") == "2"
    lines = "\n".join(stored.report_lines())
    assert "Genome annotations without sequence file: 1" in lines
    assert "Sequence level annotations (packed file, no genome file expected): 2" in lines
    for accession in ("NC_000913.3", "CP000250.1"):
        assert accession not in lines


def test_disk_report_is_replaced_by_the_next_run(tmp_path):
    db, _ = make_annotated_database(tmp_path, ["GCF_000000001.1", "GCF_000000002.1"])
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_000000001.1_genomic.fna").write_text(">a\nACGT\n")
    ProcessDirectory(db).report_disk(str(genomes))
    (genomes / "GCF_000000002.1_genomic.fna").write_text(">c\nACGT\n")
    ProcessDirectory(db).report_disk(str(genomes))
    conn = ModifyFunctions(database=db)
    assert conn.report_kinds().get("missing_on_disk") is None
    assert conn.get_meta("genomes_missing_on_disk") == "0"


def test_skipped_import_rows_are_counted_and_stay_in_the_file(tmp_path):
    db = str(tmp_path / "skip.fdb")
    read_object = ReadTaxonomy(taxonomy_file=str(tmp_path / "tax.tsv"), database=db)
    read_object.skipped = [("bad_row_1", "no taxid"), ("bad_row_2", "taxid not in database")]
    skipped_file = read_object.write_skipped()
    assert skipped_file == str(tmp_path / "tax.tsv.skipped")
    conn = ModifyFunctions(database=db)
    findings = conn.get_report("skipped_import")
    assert len(findings) == 1
    assert findings[0][0] == str(tmp_path / "tax.tsv")
    assert findings[0][1].startswith("total\t2\t")
    lines = "\n".join(conn.report_lines())
    assert "Rows skipped by the taxonomy import: 2" in lines
    assert str(skipped_file) in lines
    ## the rows themselves are only in the file, not in the database
    assert "bad_row_1" not in lines
    assert "bad_row_1\tno taxid" in open(skipped_file).read()


def test_an_older_database_still_lists_the_entries_it_stored(tmp_path):
    """Databases written by an earlier flextaxd version hold one row per finding: --report has to
    keep showing those entries instead of reading them as counts"""
    db = make_old_database(tmp_path)
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO report (kind,source,entry) VALUES ('missing_on_disk','/genomes','GCF_000000001.1')")
    conn.commit()
    lines = "\n".join(ModifyFunctions(database=db).report_lines())
    assert "Genome annotations without sequence file: 1" in lines
    assert "GCF_000000001.1" in lines
