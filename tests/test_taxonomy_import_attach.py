"""issue #162, the other face of it: the -tf import route of ReadTaxonomy.parse_taxonomy. A taxonomy
file normally uses the name "root" in its parent column, and the name then matches the root node
already in the database. Both Brucella files of the issue named their own top node instead, and the
author of the second one even wrote a self-link row "Brucella melitensis	Brucella melitensis	species"
to say "this is the top of my tree" in the way the root itself is written: the import took the row as
a level and the whole imported tree stayed off the root - 273 nodes and their links, one root that
links to nothing else. These tests keep the two things the import must now do: skip a self-link row,
and attach the parents the file never links to the database root."""

from flextaxd.modules.ReadTaxonomy import ReadTaxonomy
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions


def import_file(tmp_path, rows, name="attach.fdb"):
    taxfile = tmp_path / "tax.tsv"
    taxfile.write_text("parent\tchild\trank\n"
                       + "".join("\t".join(row) + "\n" for row in rows))
    database = str(tmp_path / name)
    ReadTaxonomy(taxonomy_file=str(taxfile), database=database).parse_taxonomy()
    return database


def test_a_file_that_names_its_own_top_node_attaches_it_to_the_root(tmp_path):
    db = import_file(tmp_path, [
        ("Genus one", "Species one", "species"),
        ("Genus one", "Species two", "species"),
    ])
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT parent FROM tree WHERE child ="
                       " (SELECT id FROM nodes WHERE name = 'Genus one')") == [(1,)]
    assert conn.validate_tree() is True


def test_a_file_that_uses_root_by_name_is_not_attached_a_second_time(tmp_path):
    """the documented shape stays exactly as it was: one link from the root, no extra attach"""
    db = import_file(tmp_path, [
        ("root", "Genus one", "no rank"),
        ("Genus one", "Species one", "species"),
    ])
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT COUNT(*) FROM tree WHERE child ="
                       " (SELECT id FROM nodes WHERE name = 'Genus one')")[0][0] == 1
    assert conn.validate_tree() is True


def test_a_self_link_row_is_skipped_and_its_node_attaches_to_the_root(tmp_path):
    """the first row of the Brucella nosub file, and the shape it was meant to say"""
    db = import_file(tmp_path, [
        ("Genus one", "Genus one", "species"),
        ("Genus one", "Species one", "species"),
    ])
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT child FROM tree WHERE parent = child AND child != 1") == []
    assert conn.select("SELECT parent FROM tree WHERE child ="
                       " (SELECT id FROM nodes WHERE name = 'Genus one')") == [(1,)]
    assert conn.validate_tree() is True


def test_a_file_with_several_top_nodes_attaches_each_of_them(tmp_path):
    """the shape the brucella_rep databases ended up with in the 0.6.0 era: three genera right under
    the root; each of them is now written by the import instead of a hand-made root row"""
    db = import_file(tmp_path, [
        ("Genus one", "Species one", "species"),
        ("Genus two", "Species two", "species"),
        ("Genus three", "Species three", "species"),
    ])
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT COUNT(*) FROM tree WHERE parent = 1 AND child != 1")[0][0] == 3
    assert conn.validate_tree() is True


def test_two_files_imported_into_one_database_join_at_the_root(tmp_path):
    """the joining people ask about: a second -tf run whose file names its top after the node the
    first file brought in. The name matches the id the database already has (#145), the two trees
    meet in that node, and the attach does not write the root link a second time - the database
    keeps exactly one self-link and one tree"""
    tax1 = tmp_path / "t1.tsv"
    tax1.write_text("parent\tchild\trank\nGenus one\tSpecies one\tspecies\n")
    tax2 = tmp_path / "t2.tsv"
    tax2.write_text("parent\tchild\trank\nGenus one\tGenus one\tspecies\nGenus one\tSpecies two\tspecies\n")
    database = str(tmp_path / "join.fdb")
    ReadTaxonomy(taxonomy_file=str(tax1), database=database).parse_taxonomy()
    ReadTaxonomy(taxonomy_file=str(tax2), database=database).parse_taxonomy()
    conn = ModifyFunctions(database=database)
    assert conn.select("SELECT COUNT(*) FROM tree WHERE parent = 1 AND child ="
                       " (SELECT id FROM nodes WHERE name = 'Genus one')")[0][0] == 1
    assert conn.select("SELECT COUNT(*) FROM tree WHERE parent = child AND child != 1")[0][0] == 0
    genus = conn.select("SELECT id FROM nodes WHERE name = 'Genus one'")[0][0]
    children = [row[0] for row in conn.select(
        "SELECT nodes.name FROM tree JOIN nodes ON nodes.id = tree.child"
        " WHERE tree.parent = ?", params=(genus,))]
    assert sorted(children) == ["Species one", "Species two"]
    assert conn.validate_tree() is True
