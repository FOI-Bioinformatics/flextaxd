"""issue #137: a CanSNPer tree adds its levels without an inherited rank.

The rows of a CanSNPer tree are strain/sequence type levels below the rank of the import, and they have no
rank of their own. The parser added the links of those levels straight through the database connection,
whose default is rank_i 1, so every one of them was labelled with whatever rank that database happened to
have added first - in a database imported from NCBI that is a real rank of the backbone. The links go
through the importer's own add_link, which takes the rank by name like every other link does.
"""

import pytest

from flextaxd.modules.ReadTaxonomyCanSNPer import ReadTaxonomyCanSNPer
from flextaxd.modules.WriteTaxonomy import WriteTaxonomy
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

"""rank rows already in the database before the import, rank_i 1 included: none of these names is ever
used by the CanSNPer importer, so a link carrying one of them was labelled by the default rank id"""
RANKS = [(1, "phylum"), (2, "genus"), (3, "species")]

"""the levels the CanSNPer file adds below the rank of the import"""
DEEP_LEVELS = ["CC17", "CC114", "ST18", "ST413", "ST308"]

TREE = """Enterococcus faecium;CC17;ST18
Enterococcus faecium;CC17;ST413
Enterococcus faecium;CC114;ST308
"""


def make_database(tmp_path):
    db = str(tmp_path / "cansnp.fdb")
    conn = ModifyFunctions(database=db)
    for rank_i, rank in RANKS:
        conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(rank_i, rank))
    conn.commit()
    return db


def make_tree(tmp_path):
    tree_file = tmp_path / "tree.txt"
    tree_file.write_text(TREE)
    return str(tree_file)


def import_taxonomy(tmp_path):
    """the root name is taken from the first field of the file, as the importer does when it is not given"""
    reader = ReadTaxonomyCanSNPer(taxonomy_file=make_tree(tmp_path), database=make_database(tmp_path),
                                  rank="species")
    reader.parse_taxonomy()
    return reader


def exported_links(reader):
    """(child, parent) -> rank field of the tree rows, read the way the dump reads them"""
    conn = ModifyFunctions(database=str(reader.database.database))
    links = {}
    for child, parent, rank in conn.select("SELECT child,parent,rank FROM tree "
                                          "LEFT JOIN rank ON rank.rank_i = tree.rank_i"):
        links[(int(child), int(parent))] = rank
    return links


def test_the_levels_of_a_cansnp_tree_do_not_inherit_a_rank_of_the_database(tmp_path):
    reader = import_taxonomy(tmp_path)
    links = exported_links(reader)

    deep = {reader.taxonomy[name] for name in DEEP_LEVELS}
    assert deep, "the tree of the test file was not added"
    for (child, parent), rank in links.items():
        if child in deep:
            """the levels of a CanSNPer file name no rank: they carry none, they do not borrow rank id 1 of
            the database (a rank of its backbone) and they are not called "no rank" either (issue #142)"""
            assert rank is None

    """the levels of the file are all there and linked to each other as the file has them"""
    children = {child for child, _ in links}
    for name in DEEP_LEVELS:
        assert reader.taxonomy[name] in children


def test_the_cansnper_levels_keep_their_shape_through_the_dump(tmp_path):
    """issue #134 kept every link of the tree; here the link that reaches the dump must not carry a rank
    name that the import never asked for"""
    reader = import_taxonomy(tmp_path)
    out = tmp_path / "dump"
    out.mkdir()
    writer = WriteTaxonomy(str(out), database=str(reader.database.database), prefix="names,nodes",
                           dbprogram="kraken2")
    writer.nodes()

    deep = {reader.taxonomy[name] for name in DEEP_LEVELS}
    with open(str(out / "nodes.dmp")) as nodes:
        ranks = {int(line.split("\t|\t")[0]): line.split("\t|\t")[2] for line in nodes if line.strip()}
    for taxid in deep:
        assert taxid in ranks
        """a level without a rank reaches the dump with an empty rank field, never the text None and never a
        rank of the backbone it was imported into (issue #134, #137, #142)"""
        assert ranks[taxid] == ""


def node_names(database):
    """every node of the database as (id, name)"""
    conn = ModifyFunctions(database=database)
    return [(int(node_id), name) for node_id, name in conn.select("SELECT id,name FROM nodes")]


def link_rows(database):
    conn = ModifyFunctions(database=database)
    return [(int(child), int(parent)) for child, parent in
            conn.select("SELECT child,parent FROM tree")]


def test_a_second_import_of_the_same_tree_adds_no_second_nodes(tmp_path):
    """issue #145: a name driven importer never looked at the nodes the database already held, so the
    second run of the same file gave every name of the tree a second taxid and the links of that run
    pointed at the copies. Nothing below the tree had two names, nothing gained two taxids."""
    database = make_database(tmp_path)
    first = ReadTaxonomyCanSNPer(taxonomy_file=make_tree(tmp_path), database=database, rank="species")
    first.parse_taxonomy()
    before_nodes = node_names(database)
    before_links = link_rows(database)

    second = ReadTaxonomyCanSNPer(taxonomy_file=make_tree(tmp_path), database=database, rank="species")
    second.parse_taxonomy()
    after_nodes = node_names(database)
    after_links = link_rows(database)

    names = [name for _, name in after_nodes]
    assert len(names) == len(set(names)), "a name of the tree got a second taxid"
    assert sorted(after_nodes) == sorted(before_nodes), "the second import added nodes"
    assert second.new_nodes == 0, "the importer counted new nodes it did not add"
    """the parent child relations of the tree did not change either: the second import hangs its levels
    on the nodes that were already there (the repeated rows a reimport writes are issue #147)"""
    assert {(child, parent) for child, parent in after_links} == \
        {(child, parent) for child, parent in before_links}, "the second import made new relations"


def test_a_tree_imported_into_a_database_that_knows_the_names_uses_them(tmp_path):
    """issue #145: the CanSNPer levels are added below the species of the database. When the species is
    already there under its own taxid, the import has to hang its levels on that node instead of adding a
    second node with the same name"""
    database = make_database(tmp_path)
    conn = ModifyFunctions(database=database)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    conn.query("INSERT INTO nodes (id,name) VALUES (2,'Enterococcus')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,2,2)")
    conn.query("INSERT INTO nodes (id,name) VALUES (3,'Enterococcus faecium')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (2,3,3)")
    conn.commit()

    reader = ReadTaxonomyCanSNPer(taxonomy_file=make_tree(tmp_path), database=database, rank="species")
    reader.parse_taxonomy()

    names = [name for _, name in node_names(database)]
    assert len(names) == len(set(names)), "a name of the tree was added a second time"
    assert sorted(names) == sorted(["root", "cellular organisms", "Enterococcus",
                                    "Enterococcus faecium"] + DEEP_LEVELS)

    """the species node of the database is the parent of the levels of the file: the name was not copied
    under a second taxid, and the tree of the import hangs on the node that was already there"""
    links = link_rows(database)
    species_taxids = [node_id for node_id, name in node_names(database)
                      if name == "Enterococcus faecium"]
    assert species_taxids == [3], "the species of the database got a second taxid: {ids}".format(ids=species_taxids)
    assert [child for child, parent in links if parent == 3], \
        "the imported levels did not hang on the species node that was already there"
