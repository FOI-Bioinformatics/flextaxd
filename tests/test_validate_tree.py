"""issue #107: validation read both whole tables into python and put every node id under root into the
text of one query, which took minutes on the 2534301 nodes of an NCBI taxdump. The checks now run as SQL
against the two tables, so the tests here are the ones that keep each of the four checks reporting what
it reported before, and one that a database holding a cycle is walked without running on for ever."""

import logging

import pytest

from flextaxd.modules.database.DatabaseConnection import ModifyFunctions, TreeError


def make_database(tmp_path, name="valid.fdb", nodes=((1, "root"), (2, "Genus one"), (3, "Species one")),
                  links=((1, 1, 1), (1, 2, 2), (2, 3, 3))):
    db = str(tmp_path / name)
    conn = ModifyFunctions(database=db)
    for rank_i, rank in ((1, "no rank"), (2, "genus"), (3, "species")):
        conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(rank_i, rank))
    for node_id, name_of_node in nodes:
        conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(node_id, name_of_node))
    for parent, child, rank_i in links:
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (?,?,?)",
                   insert_val=(parent, child, rank_i))
    conn.commit()
    return db


def test_a_tree_that_is_complete_validates(tmp_path):
    assert ModifyFunctions(make_database(tmp_path)).validate_tree() is True


def test_a_database_that_holds_nothing_validates(tmp_path):
    """the tables of a database that nothing was imported into are empty, and validation passed on the
    empty lists it read before"""
    assert ModifyFunctions(make_database(tmp_path, nodes=(), links=())).validate_tree() is True


def test_a_node_that_no_link_reaches_is_named(tmp_path):
    """a node that hangs under nothing never reaches a dump or a classifier database"""
    db = make_database(tmp_path, "orphan.fdb",
                       nodes=((1, "root"), (2, "Genus one"), (3, "Species one"), (4, "Lost node")),
                       links=((1, 1, 1), (1, 2, 2), (2, 3, 3)))
    with pytest.raises(TreeError) as error:
        ModifyFunctions(db).validate_tree()
    assert "does not match" in str(error.value)


def test_a_link_to_an_id_without_a_node_is_named(tmp_path):
    """issue #152: a taxonomy file linked a taxid that its names file never described, so the tree pointed
    at a row that does not exist; the validation found it, but only by counting the two python sets"""
    db = make_database(tmp_path, "missing.fdb",
                       nodes=((1, "root"), (2, "Genus one")),
                       links=((1, 1, 1), (1, 2, 2), (2, 999, 3)))
    with pytest.raises(TreeError) as error:
        ModifyFunctions(db).validate_tree()
    assert "annotated nodes" in str(error.value)


def test_a_link_hanging_off_a_node_that_is_not_under_root(tmp_path):
    db = make_database(tmp_path, "detached.fdb",
                       nodes=((1, "root"), (2, "Genus one"), (3, "Strain a")),
                       links=((1, 1, 1), (3, 2, 2), (2, 3, 3)))
    with pytest.raises(TreeError) as error:
        ModifyFunctions(db).validate_tree()
    assert "does not match" in str(error.value)


def test_a_cycle_under_root_is_walked_once(tmp_path):
    """a recursive WITH clause ends at the first row it has seen before: a database that holds a link back
    to an ancestor is walked once, it does not run on for ever the way a python recursion over the links
    would, and the nodes of the cycle are all under root so validation passes"""
    db = make_database(tmp_path, "cycle.fdb",
                       nodes=((1, "root"), (2, "Genus one")),
                       links=((1, 1, 1), (1, 2, 2), (2, 1, 3)))
    assert ModifyFunctions(db).validate_tree() is True


def test_the_numbers_of_the_statistics_are_the_numbers_of_the_tables(tmp_path, caplog):
    """the statistics line of a validation says what the whole tables count, which is what the length of
    the python containers said before"""
    caplog.set_level(logging.INFO)
    ModifyFunctions(make_database(tmp_path, "stats.fdb")).validate_tree()
    stats = [record.getMessage() for record in caplog.records if "Tree statistics" in record.getMessage()]
    assert len(stats) == 1
    assert "Nodes: 3" in stats[0]
    assert "Links: 3" in stats[0]
    assert "Tree: n(3), l(3)" in stats[0]
    assert "LinkNodes: 3" in stats[0]


def test_a_node_that_links_to_itself_is_named(tmp_path, caplog):
    """issue #162: the Brucella modification file started with "X X species". A self-link on a normal
    node passes every count of the other checks (the walk meets the node and stops) while its subtree
    hides from it: only the root may link to itself, and the error names the node like the other checks
    do, with the example in the log line the check writes before raising"""
    caplog.set_level(logging.INFO)
    db = make_database(tmp_path, "selflink.fdb",
                       nodes=((1, "root"), (2, "Genus one"), (3, "Species one")),
                       links=((1, 1, 1), (1, 2, 2), (2, 3, 3), (2, 2, 2)))
    with pytest.raises(TreeError) as error:
        ModifyFunctions(db).validate_tree()
    assert "linking to themselves" in str(error.value)
    assert "total: 1" in str(error.value)
    found = [record.getMessage() for record in caplog.records if "Genus one" in record.getMessage()]
    assert found, [record.getMessage() for record in caplog.records]


def test_two_root_self_links_are_named(tmp_path):
    """the root convention is exactly one self-link; brucella_rep_org.fdb holds two because the rank
    table of a pre-#138 import carried "no rank" twice and one link went to each id"""
    db = make_database(tmp_path, "rootdouble.fdb",
                       nodes=((1, "root"), (2, "Genus one")),
                       links=((1, 1, 1), (1, 2, 2), (1, 1, 2)))
    with pytest.raises(TreeError) as error:
        ModifyFunctions(db).validate_tree()
    assert "root links to itself 2 times" in str(error.value)


def test_a_doubled_rank_row_and_its_links_are_folded_when_the_database_is_opened(tmp_path):
    """the legacy shape of issue #162, built the way the old databases were: a rank table holding
    "no rank" at rank_i 1 and 2 and the root linking to itself under both ids, which the unique
    (parent,child,rank_i) index never caught because the rank_i differs. Opening the database folds
    the rank onto its lowest id, drops the copy of the root link the fold creates and lets the rest
    of the links follow, so validation can ask for exactly one root self-link and gets it"""
    import sqlite3
    db = str(tmp_path / "fold.fdb")
    old = sqlite3.connect(db)
    old.executescript(
        """CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
        CREATE TABLE tree (parent integer NOT NULL, child integer NOT NULL, rank_i integer, unique (parent,child,rank_i));
        CREATE TABLE rank (rank_i integer PRIMARY KEY, rank VARCHAR(15));
        CREATE TABLE genomes (id integer NOT NULL, genome VARCHAR(255) NOT NULL UNIQUE, reference VARCHAR(255));
        INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(2,'no rank'),(3,'genus');
        INSERT INTO nodes (id,name) VALUES (1,'root'),(2,'Genus one');
        INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1),(1,1,2),(1,2,3);""")
    old.commit()
    old.close()
    assert ModifyFunctions(database=db).validate_tree() is True
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT COUNT(*) FROM tree WHERE parent = 1 AND child = 1")[0][0] == 1
    assert sorted(row[0] for row in conn.select("SELECT rank_i FROM rank")) == [1, 3]
    assert conn.select("SELECT COUNT(*) FROM tree WHERE parent = 1 AND child = 2 AND rank_i = 3")[0][0] == 1
