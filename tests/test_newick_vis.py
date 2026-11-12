"""The visualisation options: --vis_type newick_vis and --vis_node.

issue #155: Biopython draw_ascii dies with `IndexError: list assignment index out of range` on a newick
chain of named single child nodes, which is what every taxon with exactly one child looks like, so the
ascii tree is drawn by flextaxd from the nodes it built itself.
issue #156: --vis_node only accepted a node name, although a taxid is what every other option of both
commands takes.
"""

import os

import pytest

from flextaxd.modules.NewickTree import NewickTree, VisualisationError
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

CHAIN = [(2, "Genus one", 4), (3, "Species one", 4), (4, "Strain a", 5)]


def make_chain_database(path):
    """root -> genus -> species -> strain: three levels that each have exactly one child"""
    conn = ModifyFunctions(database=str(path))
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(4,'species'),(5,'strain')")
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    for taxid, name, rank in CHAIN:
        conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(taxid, name))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (?,?,?)",
                   insert_val=(taxid - 1, taxid, rank))
    conn.commit()
    return str(path)


def test_ascii_tree_prints_a_chain_of_single_children(tmp_path, capsys, monkeypatch):
    """the shape draw_ascii cannot print: `(((Strain a)Species one)Genus one)root;`"""
    monkeypatch.chdir(tmp_path)
    database = make_chain_database(tmp_path / "vis.fdb")
    tree = NewickTree(database=database, taxid="root", maxdepth=0)
    assert "(((" in tree.newickTree
    assert tree.print("newick_vis") is True
    printed = capsys.readouterr().out
    for name in ("root", "Genus one", "Species one", "Strain a"):
        assert name in printed
    assert "`-- " in printed and "|-- " not in printed
    assert os.path.exists(str(tmp_path / "flextaxd.vis.png"))


def test_ascii_tree_prints_a_node_with_two_children_once_each(tmp_path, capsys, monkeypatch):
    """a node with two children prints both branches, each once"""
    monkeypatch.chdir(tmp_path)
    database = make_chain_database(tmp_path / "two.fdb")
    conn = ModifyFunctions(database=database)
    conn.query("INSERT INTO nodes (id,name) VALUES (5,'Strain b')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (3,5,5)")
    conn.commit()
    tree = NewickTree(database=database, taxid="root", maxdepth=0)
    tree.print("newick_vis")
    printed = capsys.readouterr().out
    assert printed.count("Strain a") == 1 and printed.count("Strain b") == 1
    assert "|-- Strain a" in printed and "`-- Strain b" in printed


def test_vis_node_accepts_a_taxid(tmp_path):
    database = make_chain_database(tmp_path / "taxid.fdb")
    assert NewickTree(database=database, taxid=2, maxdepth=0).taxid == 2
    assert NewickTree(database=database, taxid="2", maxdepth=0).taxid == 2


def test_unknown_taxid_is_named(tmp_path):
    """a taxid that is not in the database gets the error naming it, not the name lookup failure"""
    database = make_chain_database(tmp_path / "unknown.fdb")
    with pytest.raises(VisualisationError):
        NewickTree(database=database, taxid=999999, maxdepth=0)
