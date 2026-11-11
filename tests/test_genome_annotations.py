import os
import sqlite3

import pytest

from flextaxd.modules import ModifyTree as modify_module
from flextaxd.modules.ReadTaxonomy import ReadTaxonomy
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions


def make_database(path):
    """Create a database with a root and two species nodes, without genome annotations"""
    db = str(path)
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(4,'species')")
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    conn.query("INSERT INTO nodes (id,name) VALUES (2,'Escherichia coli')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,2,4)")
    conn.query("INSERT INTO nodes (id,name) VALUES (3,'Salmonella enterica')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,3,4)")
    conn.commit()
    return db


def genome_rows(db):
    return {genome: (node_i, reference) for node_i, genome, reference in
            ModifyFunctions(database=db).select("SELECT id,genome,reference FROM genomes")}


def write_file(path, rows):
    file = str(path)
    with open(file, "w") as out:
        out.write("".join("\t".join(row) + "\n" for row in rows))
    return file


def annotate(db, annotation_file, reference=False):
    """Run the annotation route a lone --genomeid2taxid gives (custom_taxonomy_databases.py:372)"""
    modify_object = modify_module.ModifyTree(database=db, update_genomes=True)
    modify_object.update_annotations(genomeid2taxid=annotation_file, reference=reference)


def test_annotation_file_reference_column_is_stored(tmp_path):
    """Regression test: the third column of the annotation file was unpacked and dropped, so the
    modify route left every genome without its refseq/genbank label while the create route stored it,
    and the genome annotation dump (genome, name, reference) could not be read back losslessly"""
    db = make_database(tmp_path / "ann.fdb")
    annotation_file = write_file(tmp_path / "ann.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "refseq"),
        ("GCA_000000002.1", "Escherichia coli", "genbank"),
        ("GCF_000000003.1", "Salmonella enterica", "refseq"),
    ])
    annotate(db, annotation_file)
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "refseq")
    assert rows["GCA_000000002.1"] == (2, "genbank")
    assert rows["GCF_000000003.1"] == (3, "refseq")


def test_reference_argument_overrides_the_column_of_the_file(tmp_path):
    """A reference given as argument labels every genome, as in the create route where the column of
    the file is then only read as an override"""
    db = make_database(tmp_path / "override.fdb")
    annotation_file = write_file(tmp_path / "override.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "genbank"),
        ("GCF_000000002.1", "Salmonella enterica", "genbank"),
    ])
    annotate(db, annotation_file, reference="nt")
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "nt")
    assert rows["GCF_000000002.1"] == (3, "nt")


def test_rows_without_reference_are_not_labelled_from_the_previous_row(tmp_path):
    """A two column row keeps the reference empty, it must not inherit the value of the row before
    it (the same walking state that mislabelled genomes in issue #113)"""
    db = make_database(tmp_path / "plain.fdb")
    annotation_file = write_file(tmp_path / "plain.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "genbank"),
        ("GCF_000000002.1", "Escherichia coli"),
    ])
    annotate(db, annotation_file)
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "genbank")
    assert rows["GCF_000000002.1"] == (2, None)


def read_annotate(db, annotation_file, reference=False):
    """Run the annotation route a --genomeid2taxid during a build gives (ReadTaxonomy)"""
    taxonomy_file = os.path.join(os.path.dirname(str(db)), "taxonomy.tsv")
    read_object = ReadTaxonomy(taxonomy_file=taxonomy_file, database=db)
    read_object.parse_genomeid2taxid(annotation_file, reference=reference)
    return read_object.skipped


def test_create_route_reference_does_not_carry_over(tmp_path):
    """Regression test: the create route reused its reference argument as state between rows, so a
    two column row was labelled with the reference of the last row that had a third column"""
    db = make_database(tmp_path / "create.fdb")
    annotation_file = write_file(tmp_path / "create.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "genbank"),
        ("GCF_000000002.1", "Escherichia coli"),
        ("GCF_000000003.1", "Salmonella enterica", "refseq"),
    ])
    assert read_annotate(db, annotation_file) == []
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "genbank")
    assert rows["GCF_000000002.1"] == (2, None)
    assert rows["GCF_000000003.1"] == (3, "refseq")


def test_create_route_extra_column_does_not_abort_the_import(tmp_path):
    """Regression test: a row with four columns passed the two column unpack only to die in the
    three column one, an uncaught ValueError that ended the whole import"""
    db = make_database(tmp_path / "extra.fdb")
    annotation_file = write_file(tmp_path / "extra.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "genbank", "extra"),
        ("GCF_000000002.1", "Salmonella enterica", "refseq"),
    ])
    assert read_annotate(db, annotation_file) == []
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "genbank")
    assert rows["GCF_000000002.1"] == (3, "refseq")


def test_create_route_row_without_two_columns_is_reported_not_fatal(tmp_path):
    """A row that cannot give a genome and a node is listed as skipped instead of crashing"""
    db = make_database(tmp_path / "broken.fdb")
    annotation_file = write_file(tmp_path / "broken.tsv", [
        ("only_one_column",),
        ("GCF_000000001.1", "Escherichia coli", "genbank"),
    ])
    skipped = read_annotate(db, annotation_file)
    assert [genome for genome, reason in skipped] == ["only_one_column"]
    assert genome_rows(db) == {"GCF_000000001.1": (2, "genbank")}


def test_create_route_reference_argument_overrides_the_column(tmp_path):
    """The reference given as argument labels every row, the column of the file is then ignored"""
    db = make_database(tmp_path / "createref.fdb")
    annotation_file = write_file(tmp_path / "createref.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "genbank"),
        ("GCF_000000002.1", "Salmonella enterica"),
    ])
    assert read_annotate(db, annotation_file, reference="nt") == []
    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "nt")
    assert rows["GCF_000000002.1"] == (3, "nt")


def test_a_taxid_the_database_has_no_node_for_is_skipped(tmp_path, caplog):
    """issue #73: a row whose node is a number was taken as a taxid and never asked of the database, so a
    file naming a taxid the database does not have pointed the genome at a node that does not exist. The
    genome then vanished from the dump and from the kraken2 build, and the run said nothing about it,
    while the same situation under a node NAME was already listed in the skipped file"""
    db = make_database(tmp_path / "taxid.fdb")
    annotation_file = write_file(tmp_path / "taxid.tsv", [
        ("GCF_000000001.1", "2", "refseq"),
        ("GCF_000000002.1", "999999", "refseq"),
    ])
    annotate(db, annotation_file)

    rows = genome_rows(db)
    assert rows["GCF_000000001.1"] == (2, "refseq"), "a taxid that is in the database stopped working"
    assert "GCF_000000002.1" not in rows, "a genome was pointed at a node that does not exist"

    with open(annotation_file + ".skipped") as skipped:
        listed = skipped.read()
    assert "GCF_000000002.1" in listed and "999999" in listed
    warnings = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    assert any("999999" in warning for warning in warnings), warnings
    assert ModifyFunctions(database=db).get_report("skipped_import", annotation_file) == [
        (annotation_file, "total\t1\t" + annotation_file + ".skipped")]


def test_an_unknown_taxid_is_skipped_by_the_batched_path_as_well(tmp_path):
    """the check is not tied to the size of the file: a large file takes the batched write and still loses
    the rows whose taxid has no node, before anything is written"""
    db = make_database(tmp_path / "bigtaxid.fdb")
    annotation_file = write_file(tmp_path / "bigtaxid.tsv", [
        ("GCF_{:09d}.1".format(i),
         "999999" if i % 10 == 0 else ("2" if i % 2 == 0 else "Salmonella enterica"),
         "refseq")
        for i in range(150)
    ])
    annotate(db, annotation_file)

    rows = genome_rows(db)
    assert len(rows) == 135, "the rows naming an unknown taxid were written anyway: {n}".format(n=len(rows))
    for i in range(150):
        genome = "GCF_{:09d}.1".format(i)
        if i % 10 == 0:
            assert genome not in rows
        else:
            assert rows[genome][0] == (2 if i % 2 == 0 else 3)
    with open(annotation_file + ".skipped") as skipped:
        assert len(skipped.read().strip().splitlines()) == 15


def test_a_large_annotation_file_is_written_in_batches(tmp_path, monkeypatch):
    """issue #126: the file driven route had no batched path, every row cost one or two single UPDATE
    statements, about 800000 statements for the 394932 genomes of the GTDB case, where the merge route
    does it in one (multi_update, issue #116). Above the 100 row threshold of the merge route the file
    route now batches too, and keeps the same numbers of added and updated annotations."""
    db = make_database(tmp_path / "many.fdb")
    annotation_file = write_file(tmp_path / "many.tsv", [
        ("GCF_{:09d}.1".format(i),
         "2" if i % 2 == 0 else "Salmonella enterica",
         "refseq" if i % 3 else "genbank")
        for i in range(250)
    ])

    single = []
    original_update_genome = ModifyFunctions.update_genome
    def spy_update_genome(self, data):
        single.append(data["where"])
        return original_update_genome(self, data)
    monkeypatch.setattr(ModifyFunctions, "update_genome", spy_update_genome)

    batches = []
    original_multi_update = ModifyFunctions.multi_update
    def spy_multi_update(self, data, table):
        batches.append(len(data["data"]))
        return original_multi_update(self, data, table)
    monkeypatch.setattr(ModifyFunctions, "multi_update", spy_multi_update)

    annotate(db, annotation_file)

    assert single == [], "{n} single UPDATE statements were still run".format(n=len(single))
    assert batches and sum(batches) == 250
    rows = genome_rows(db)
    assert len(rows) == 250
    for i in range(250):
        genome = "GCF_{:09d}.1".format(i)
        assert rows[genome] == (2 if i % 2 == 0 else 3, "refseq" if i % 3 else "genbank")


def test_a_small_annotation_file_keeps_the_row_by_row_path(tmp_path, monkeypatch):
    """the numbers of the row by row path (a genome that is not in the table yet is added, the others are
    updated) are kept for the small files that are annotated one at a time"""
    db = make_database(tmp_path / "few.fdb")
    seed = ModifyFunctions(database=db)
    seed.query("INSERT INTO genomes (id,genome) VALUES (2,'GCF_000000001.1')")
    seed.commit()
    annotation_file = write_file(tmp_path / "few.tsv", [
        ("GCF_000000001.1", "Escherichia coli", "refseq"),
        ("GCF_000000002.1", "Salmonella enterica", "genbank"),
    ])

    batches = []
    monkeypatch.setattr(ModifyFunctions, "multi_update",
                        lambda self, data, table: batches.append(len(data["data"])) or True)
    annotate(db, annotation_file)

    assert batches == [], "a file of two rows was batched"
    assert genome_rows(db) == {"GCF_000000001.1": (2, "refseq"), "GCF_000000002.1": (3, "genbank")}
