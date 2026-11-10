import logging
import os
import sqlite3

import pytest

from flextaxd.modules import ModifyTree as modify_module
from flextaxd.modules.database.DatabaseConnection import DatabaseFunctions, ModifyFunctions


def make_database(path, taxon, taxid, genomes):
    """Create a database with a root, one taxon node below it and one annotated genome per entry,
    every annotation carrying its reference label"""
    db = str(path)
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(4,'species')")
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(taxid, taxon))
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,{},4)".format(taxid))
    for genome, reference in genomes:
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, genome, reference))
    conn.commit()
    return db


def merge(target, source):
    """Run the merge a --mod_database gives: the branch of the source database is added to the
    receiving database and its genome annotations are transfered"""
    modify_object = modify_module.ModifyTree(database=target, mod_database=source, parent="root")
    modify_object.update_database()


def genome_rows(db):
    return {genome: (node_i, reference) for node_i, genome, reference in
            ModifyFunctions(database=db).select("SELECT id,genome,reference FROM genomes")}


def test_merge_transfers_the_reference_of_genomes(tmp_path):
    """Regression test: the merge used to write only the id column, so transferred genomes lost
    their refseq/genbank label (GTDB_NCBI_ftul_merge.db: 394932 rows with reference NULL)"""
    target = make_database(tmp_path / "target.fdb", "Salmonella enterica", 2,
                           [("GCF_000000001.1", "refseq")])
    source = make_database(tmp_path / "source.fdb", "Escherichia coli", 412,
                           [("GCF_000000002.1", "refseq"), ("GCA_000000003.1", "genbank")])
    merge(target, source)
    rows = genome_rows(target)
    assert rows["GCF_000000001.1"][1] == "refseq"
    assert rows["GCF_000000002.1"][1] == "refseq"
    assert rows["GCA_000000003.1"][1] == "genbank"
    ## the transferred genomes point at the node the receiving database made for the taxon
    new_node = ModifyFunctions(database=target).get_node("Escherichia coli")[0]
    assert new_node != 412
    assert rows["GCF_000000002.1"][0] == new_node
    assert rows["GCA_000000003.1"][0] == new_node


def test_batched_merge_transfers_the_reference_of_genomes(tmp_path):
    """The same for a merge big enough to take the batched update path (100 or more genomes)"""
    genomes = [("GCF_{:09d}.1".format(i), "genbank" if i % 2 else "refseq") for i in range(2, 150)]
    target = make_database(tmp_path / "target.fdb", "Salmonella enterica", 2,
                           [("GCF_000000001.1", "refseq")])
    source = make_database(tmp_path / "source.fdb", "Escherichia coli", 412, genomes)
    merge(target, source)
    rows = genome_rows(target)
    assert len(rows) == len(genomes) + 1
    assert [reference for _, reference in rows.values()].count(None) == 0
    for genome, reference in genomes:
        assert rows[genome][1] == reference


def test_multi_update_keeps_a_label_when_the_incoming_value_is_null(tmp_path):
    """A source database without a reference for one of its genomes does not erase the label the
    receiving database already stored for that genome"""
    db = str(tmp_path / "null.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (7,'GCF_000000001.1','refseq')")
    conn.commit()
    conn.multi_update({"set_column": "id", "where_column": "genome",
                       "columns": ["id", "genome", "reference"],
                       "data": ["(9,'GCF_000000001.1',NULL)"]}, "genomes")
    conn.commit()
    assert genome_rows(db) == {"GCF_000000001.1": (9, "refseq")}


def test_multi_update_updates_existing_rows(tmp_path):
    """The batched update used to assign the existing row to itself, so an existing row never got
    the incoming values: a genome that is already annotated moves to the new node and keeps the
    reference that travels with it"""
    db = str(tmp_path / "batch.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (7,'GCF_000000001.1','refseq')")
    conn.commit()
    conn.multi_update({"set_column": "id", "where_column": "genome",
                       "columns": ["id", "genome", "reference"],
                       "data": ["(9,'GCF_000000001.1','genbank')", "(9,'GCF_000000002.1','genbank')"]},
                      "genomes")
    conn.commit()
    assert genome_rows(db) == {"GCF_000000001.1": (9, "genbank"), "GCF_000000002.1": (9, "genbank")}


def test_multi_update_reports_a_failed_statement(tmp_path):
    """Regression test: the batched update returned True whatever the query did, so a merge whose
    transfer failed was reported as complete"""
    db = str(tmp_path / "failed.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (7,'GCF_000000001.1','refseq')")
    conn.commit()
    assert conn.multi_update({"set_column": "id", "where_column": "genome",
                              "columns": ["id", "genome", "reference"],
                              "data": ["(9,'GCF_000000001.1','genbank')"]}, "not_a_table") is False
    assert genome_rows(db) == {"GCF_000000001.1": (7, "refseq")}


def test_batched_merge_reports_the_numbers_transfered(tmp_path, caplog):
    """The batched path reports the same numbers as the per-genome path instead of saying that no
    statistics were calculated"""
    caplog.set_level(logging.INFO)
    genomes = [("GCF_{:09d}.1".format(i), "refseq") for i in range(2, 150)]
    target = make_database(tmp_path / "count.fdb", "Salmonella enterica", 2,
                           [("GCF_000000001.1", "refseq")])
    source = make_database(tmp_path / "countsrc.fdb", "Escherichia coli", 412, genomes)
    merge(target, source)
    assert "148 added and 0 genome annotations were transfered in one batch" in caplog.text


def test_batched_merge_reports_a_failed_transfer(tmp_path, caplog, monkeypatch):
    """A batched transfer that the database refused is named as an error and the receiving database
    keeps the annotations it had"""
    caplog.set_level(logging.INFO)
    monkeypatch.setattr(DatabaseFunctions, "multi_update", lambda self, data, table: False)
    genomes = [("GCF_{:09d}.1".format(i), "refseq") for i in range(2, 150)]
    target = make_database(tmp_path / "failmerge.fdb", "Salmonella enterica", 2,
                           [("GCF_000000001.1", "refseq")])
    source = make_database(tmp_path / "failsrc.fdb", "Escherichia coli", 412, genomes)
    merge(target, source)
    assert "148 genome annotations could not be transfered" in caplog.text
    assert genome_rows(target) == {"GCF_000000001.1": (2, "refseq")}


def test_merge_lists_genomes_it_could_not_transfer_in_a_file(tmp_path):
    """A genome whose node does not come over with the merged branch is listed in a file next to the
    receiving database, the database itself keeps only the number"""
    target = make_database(tmp_path / "skipmerge.fdb", "Salmonella enterica", 2,
                           [("GCF_000000001.1", "refseq")])
    source = make_database(tmp_path / "skipsrc.fdb", "Escherichia coli", 412,
                           [("GCF_000000002.1", "refseq")])
    conn = ModifyFunctions(database=source)
    conn.query("INSERT INTO nodes (id,name) VALUES (900,'Clostridium botulinum')")
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (900,'GCF_000000009.1','genbank')")
    conn.commit()
    merge(target, source)
    notadded_file = target + ".notadded"
    assert os.path.exists(notadded_file)
    listed = open(notadded_file).read()
    assert "GCF_000000009.1" in listed and "Clostridium botulinum" in listed
    findings = ModifyFunctions(database=target).get_report("skipped_import")
    assert len(findings) == 1
    assert findings[0][1].startswith("total\t1\t")
    entries = " ".join(entry for _, entry in findings)
    assert "GCF_000000009.1" not in entries
