import os
import sqlite3

import pytest

from flextaxd.modules import ModifyTree as modify_module
from flextaxd.modules.ProcessDirectory import ProcessDirectory
from flextaxd.modules.database import DatabaseConnection
from flextaxd.modules.database.DatabaseConnection import DatabaseFunctions, ModifyFunctions


def make_database(tmp_path, n_genomes=0, n_nodes=10):
    """Create an empty FlexTaxD database and fill it with nodes and genome annotations"""
    db = str(tmp_path / "test.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    ## FlexTaxD databases carry a self link on the root node, validate_tree expects it
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    node_ids = [i for i in range(2, 2 + n_nodes)]
    conn.query("INSERT INTO nodes (id,name) VALUES ({},'root_child')".format(node_ids[0]))
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,{},1)".format(node_ids[0]))
    for node in node_ids[1:]:
        conn.query("INSERT INTO nodes (id,name) VALUES ({},'leaf_{}')".format(node, node))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES ({},{},1)".format(node_ids[0], node))
    genomes = ["GENOME_{:07d}".format(i) for i in range(n_genomes)]
    conn.cursor.executemany(
        "INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
        [(node_ids[i % len(node_ids)], genome, "refseq") for i, genome in enumerate(genomes)],
    )
    conn.commit()
    return db, genomes, node_ids


def row_count(db, table):
    return sqlite3.connect(db).execute("SELECT count(*) FROM {}".format(table)).fetchone()[0]


def test_delete_genomes_small_list(tmp_path):
    db, genomes, _ = make_database(tmp_path, n_genomes=50)
    conn = ModifyFunctions(database=db)
    assert conn.delete_genomes([], genomes=genomes[:10], match_genome_only=True) is True
    assert row_count(db, "genomes") == 40


def test_delete_genomes_large_list(tmp_path):
    """Regression test: purging more genomes than sqlite accepts as bind variables (250000 on
    debian/ubuntu builds, 32766 on upstream builds) must still remove every genome"""
    db, genomes, _ = make_database(tmp_path, n_genomes=300000, n_nodes=100)
    conn = ModifyFunctions(database=db)
    assert conn.delete_genomes([], genomes=set(genomes), match_genome_only=True) is True
    assert row_count(db, "genomes") == 0


def test_delete_genomes_large_list_partial(tmp_path):
    db, genomes, _ = make_database(tmp_path, n_genomes=300000, n_nodes=100)
    conn = ModifyFunctions(database=db)
    keep = set(genomes[:44164])
    assert conn.delete_genomes([], genomes=set(genomes[44164:]), match_genome_only=True) is True
    assert row_count(db, "genomes") == len(keep)


def test_delete_genomes_nodes_and_genomes(tmp_path):
    db, genomes, node_ids = make_database(tmp_path, n_genomes=50000, n_nodes=100)
    conn = ModifyFunctions(database=db)
    nodes = node_ids[:5]
    owned = [g for (g,) in conn.query(
        "SELECT genome FROM genomes WHERE id in({})".format(",".join(map(str, nodes)))).fetchall()]
    assert len(owned) > DatabaseConnection.MAX_BIND_PARAMS
    assert conn.delete_genomes(nodes, genomes=owned) is True
    assert row_count(db, "genomes") == 50000 - len(owned)


def test_delete_genomes_reports_failure(tmp_path):
    db, genomes, _ = make_database(tmp_path, n_genomes=5000)
    conn = ModifyFunctions(database=db)
    assert conn.delete_genomes([], genomes=genomes, table="not_a_table") is False
    assert row_count(db, "genomes") == 5000


def test_fast_delete_links_exact(tmp_path):
    """Regression test: links are removed as exact (parent,child,rank) combinations, a link that
    combines a removed parent with a kept child must survive"""
    db = str(tmp_path / "links.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO nodes (id,name) VALUES (2,'p1')")
    conn.query("INSERT INTO nodes (id,name) VALUES (3,'p2')")
    conn.query("INSERT INTO nodes (id,name) VALUES (4,'c1')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,2,1)")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,3,1)")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (2,4,4)")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (3,4,4)")
    conn.commit()
    assert conn.fast_delete_links([(2, 4, 4), (1, 2, 1)]) is True
    assert conn.query("SELECT parent,child FROM tree").fetchall() == [(1, 3), (3, 4)]


def test_fast_delete_links_large_list(tmp_path):
    db, _, _ = make_database(tmp_path, n_genomes=0, n_nodes=30000)
    conn = ModifyFunctions(database=db)
    links = conn.query("SELECT parent,child,rank_i FROM tree WHERE parent != 1").fetchall()
    assert len(links) > DatabaseConnection.MAX_BIND_PARAMS
    assert conn.fast_delete_links(links) is True
    ## the root self link and the link to the first node hang off the root and are kept
    assert row_count(db, "tree") == 2


def test_delete_nodes_large_list(tmp_path):
    db, _, node_ids = make_database(tmp_path, n_genomes=0, n_nodes=30000)
    conn = ModifyFunctions(database=db)
    assert conn.delete_nodes(node_ids) is True
    assert row_count(db, "nodes") == 1


def test_delete_nodes_reports_failure(tmp_path):
    db, _, node_ids = make_database(tmp_path, n_genomes=0, n_nodes=3000)
    conn = ModifyFunctions(database=db)
    assert conn.delete_nodes(node_ids, table="not_a_table") is False


def test_purge_aborts_cleaning_when_delete_fails(tmp_path, monkeypatch):
    db, genomes, _ = make_database(tmp_path, n_genomes=2000)
    cleaned = []
    monkeypatch.setattr(DatabaseFunctions, "delete_genomes", lambda self, *a, **kw: False)
    monkeypatch.setattr(
        modify_module.ModifyTree, "clean_database", lambda self, ncbi=False: cleaned.append(True))
    modify_object = modify_module.ModifyTree(database=db, purge_database=True)
    modify_object.purge_database(set(genomes))
    assert cleaned == []
    assert row_count(db, "genomes") == 2000


def make_branch(conn, nodes, parent=1):
    """Insert a chain of nodes below a parent node, creating the root node when missing"""
    if parent == 1 and not conn.query("SELECT id FROM nodes WHERE id = 1").fetchone():
        conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    for node in nodes:
        conn.query("INSERT INTO nodes (id,name) VALUES ({},'n{}')".format(node, node))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES ({},{},4)".format(parent, node))
        parent = node


def row_ids(db, table, column="id"):
    return set(row[0] for row in sqlite3.connect(db).execute("SELECT {} FROM {}".format(column, table)))


def test_ancestor_and_descendant_nodes(tmp_path):
    db = str(tmp_path / "walk.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5])
    make_branch(conn, [10, 11])
    conn.commit()
    assert set(conn.get_ancestor_nodes([5])) == {1, 2, 3, 4, 5}
    assert set(conn.get_descendant_nodes([2], max_depth=2)) == {3, 4}
    assert set(conn.get_descendant_nodes([1], max_depth=1)) == {2, 10}
    assert conn.get_genome_nodes() == []


def test_delete_rows_not_in(tmp_path):
    db, _, _ = make_database(tmp_path, n_genomes=0, n_nodes=3000)
    conn = ModifyFunctions(database=db)
    temp_table = conn._temp_value_table("ftx_test_keep", [1, 2, 3])
    assert temp_table
    assert conn.delete_rows_not_in("nodes", "id", temp_table) == 2998
    assert row_ids(db, "nodes") == {1, 2, 3}


def test_clean_database_keeps_annotated_and_top_levels(tmp_path):
    """Annotated branches and their parents are kept, the top three levels of every branch are
    kept so that new taxonomy can be merged in later, everything else is removed"""
    db = str(tmp_path / "clean.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5])            ## annotated, five levels deep
    make_branch(conn, [10, 11, 12, 13, 14])    ## never annotated
    make_branch(conn, [20, 21, 22, 23, 24])    ## annotated at the deepest node only
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (5,'g1','refseq')")
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (24,'g2','refseq')")
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, clean_database=True)
    modify_object.clean_database(ncbi=True)
    assert row_ids(db, "nodes") == {1, 2, 3, 4, 5, 10, 11, 12, 20, 21, 22, 23, 24}
    assert row_ids(db, "tree", "child") == {1, 2, 3, 4, 5, 10, 11, 12, 20, 21, 22, 23, 24}
    assert row_count(db, "genomes") == 2


def test_clean_database_non_ncbi_keeps_annotated_branch(tmp_path):
    """Without NCBI mode the top levels are not kept, only annotated nodes and their parents"""
    db = str(tmp_path / "clean2.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5, 6, 7])
    make_branch(conn, [10, 11])
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (7,'g1','refseq')")
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, clean_database=True)
    modify_object.clean_database(ncbi=False)
    assert row_ids(db, "nodes") == {1, 2, 3, 4, 5, 6, 7}
    assert row_count(db, "genomes") == 1


def test_clean_database_keeps_root_without_annotated_parents(tmp_path):
    """Regression test: the root node must survive an NCBI clean also when it is not a parent of an
    annotated node, otherwise the links hanging off the root point at a deleted node and the
    validation of the cleaned database fails"""
    db = str(tmp_path / "clean3.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3])            ## below the root, never annotated
    conn.query("INSERT INTO nodes (id,name) VALUES (10,'disconnected')")
    conn.query("INSERT INTO nodes (id,name) VALUES (11,'annotated')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (10,11,4)")
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (11,'g1','refseq')")
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, clean_database=True)
    modify_object.clean_database(ncbi=True)
    assert row_ids(db, "nodes") == {1, 2, 3, 10, 11}
    assert row_ids(db, "tree", "parent") == {1, 2, 10}


def test_clean_database_without_annotations_keeps_top_levels(tmp_path):
    """An NCBI clean of an unannotated database is allowed, it is used to slim a database before
    adding the taxonomy that is wanted: the root and its three levels survive and the size of the
    result is warned about"""
    db = str(tmp_path / "clean4.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5, 6])
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, clean_database=True)
    modify_object.clean_database(ncbi=True)
    assert row_ids(db, "nodes") == {1, 2, 3, 4}
    assert row_count(db, "tree") == 4


def test_purge_removes_missing_genomes_from_shared_nodes(tmp_path):
    """The removed --purge_database_force option needs no replacement: a purge matches genomes on
    their name, so the annotation of a missing genome is also removed from a node that keeps a
    genome which is present, the GTDB case where one species node carries the representative and
    all other genomes of that species"""
    db = str(tmp_path / "purge1.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5, 6])
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (6,'g_present','refseq')")
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (6,'g_missing','genbank')")
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, purge_database=True)
    modify_object.purge_database({"g_missing"})
    assert [g for (g,) in conn.query("SELECT genome FROM genomes").fetchall()] == ["g_present"]
    assert row_ids(db, "nodes") == {1, 2, 3, 4, 5, 6}


def test_purge_removes_nodes_without_annotated_genome_below_them(tmp_path):
    """A purge cleans the database afterwards: a branch whose only genome is missing loses its nodes
    below the top levels that the clean keeps, the branch with a genome that is present survives"""
    db = str(tmp_path / "purge2.fdb")
    conn = ModifyFunctions(database=db)
    make_branch(conn, [2, 3, 4, 5, 6])         ## only genome of this branch is missing
    make_branch(conn, [10, 11, 12, 13, 14])    ## genome of this branch is present
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (6,'g_missing','genbank')")
    conn.query("INSERT INTO genomes (id,genome,reference) VALUES (14,'g_present','refseq')")
    conn.commit()
    modify_object = modify_module.ModifyTree(database=db, purge_database=True)
    modify_object.purge_database({"g_missing"})
    assert row_count(db, "genomes") == 1
    assert row_ids(db, "nodes") == {1, 2, 3, 4, 10, 11, 12, 13, 14}
    assert row_ids(db, "tree", "child") == {1, 2, 3, 4, 10, 11, 12, 13, 14}


def make_gtdb_shape_database(path, n_species=200, genomes_per_species=20):
    """One node per species carrying its representative and the other genomes of that species, the
    shape of the full GTDB taxonomy with only the representative dataset on disk"""
    db = str(path)
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(4,'species')")
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    for n in range(n_species):
        node = 2 + n
        conn.query("INSERT INTO nodes (id,name) VALUES ({},'species_{:04d}')".format(node, n))
    rows = []
    for n in range(n_species):
        for g in range(genomes_per_species):
            rows.append((2 + n, "GCF_{:09d}.1".format(n * 100 + g), "refseq" if g == 0 else "genbank"))
    conn.cursor.executemany("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)", rows)
    conn.commit()
    representatives = [rows[n * genomes_per_species][1] for n in range(n_species)]
    return db, [row[1] for row in rows], representatives


def test_purge_does_not_grow_the_database_and_leaves_no_stale_report(tmp_path):
    """Regression test for the two effects that made --purge_database enlarge the database: the scan
    of the genome folder stored one report row per missing genome, and the pages the purge freed
    stayed in the file because the clean had no nodes to remove (every species keeps its
    representative) and therefore skipped its vacuum"""
    db, accessions, representatives = make_gtdb_shape_database(
        tmp_path / "gtdb.fdb", n_species=300, genomes_per_species=30)
    folder = tmp_path / "genomes"
    folder.mkdir()
    for accession in representatives:
        (folder / (accession + "_ASMv1_genomic.fna")).write_text(">x\nACGT\n")
    size_before = os.path.getsize(db)

    files, missing = ProcessDirectory(db).process_folder(str(folder))
    missing_genomes = set(entry["genome_id"] for entry in missing)
    assert len(missing_genomes) == len(accessions) - len(representatives)
    ## the scan keeps the number, not one row per missing genome
    assert row_count(db, "report") == 1
    assert ModifyFunctions(database=db).get_meta("genomes_missing_on_disk") == str(len(missing_genomes))

    modify_object = modify_module.ModifyTree(database=db, purge_database=True)
    modify_object.purge_database(missing_genomes, genome_path=str(folder))

    assert row_count(db, "genomes") == len(representatives)
    conn = ModifyFunctions(database=db)
    assert conn.get_meta("genomes_missing_on_disk") == "0"
    assert conn.get_report("missing_on_disk") == []
    stored = " ".join(entry for entry, in conn.select("SELECT entry FROM report"))
    removed = set(accessions) - set(representatives)
    assert not [accession for accession in removed if accession in stored]
    ## the vacuum gives the freed pages back, so the command cannot leave the file bigger behind
    assert conn.query("PRAGMA freelist_count").fetchone()[0] == 0
    assert os.path.getsize(db) <= size_before

