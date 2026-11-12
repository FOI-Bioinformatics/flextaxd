"""issue #127: the modification routes read their input files with plain open()

The three files the modify route takes from the command line are the same kind of input as the files the
build route reads, and some of them are routinely compressed: the genome annotation file
(`--genomeid2taxid`), the modification file (`-mf`) and the name change file (`--update_names`). The build
route opens them through the gz sniffing opener shared by the builders and the taxonomy readers (issue
#117), so a compressed file works there and reads as binary garbage here: the row loop of a modification
file never finds a tab, and the annotation file is not parsed at all.

The second thing these tests pin is that the modification file is read row by row: the route is the file
version of the merge route, which walks every link of the source database.
"""

import gzip

from flextaxd.modules import ModifyTree as modify_module
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

ANNOTATION_ROWS = [
    ("GCF_000000001.1", "Salmonella enterica", "refseq"),
    ("GCA_000000002.1", "Shigella dysenteriae", "genbank"),
]
MODIFICATION_ROWS = [
    ("parent", "child", "level"),
    ("root", "Enterobacteriaceae", "family"),
    ("Enterobacteriaceae", "Escherichia coli", "species"),
]


def make_database(path):
    """a root with one species node below it, no genome annotations"""
    db = str(path)
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (1,'no rank'),(4,'species'),(5,'family')")
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    conn.query("INSERT INTO nodes (id,name) VALUES (2,'Salmonella enterica')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,2,4)")
    conn.query("INSERT INTO nodes (id,name) VALUES (3,'Shigella dysenteriae')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,3,4)")
    conn.commit()
    return db


def write(path, rows, compressed=False):
    """the same rows, plain or gzipped"""
    name = str(path)
    body = "".join("\t".join(row) + "\n" for row in rows)
    handle = gzip.open(name, "wt") if compressed else open(name, "w")
    with handle:
        handle.write(body)
    return name


def links(db):
    return {(parent, child) for parent, child in
            ModifyFunctions(database=db).select("SELECT parent,child FROM tree") if parent != child}


def node_names(db):
    return {name for _, name in ModifyFunctions(database=db).select("SELECT id,name FROM nodes")}


def genome_rows(db):
    return {genome: (node_i, reference) for node_i, genome, reference in
            ModifyFunctions(database=db).select("SELECT id,genome,reference FROM genomes")}


def test_a_compressed_annotation_file_annotates_like_a_plain_one(tmp_path):
    """the modify route of --genomeid2taxid must read the .gz file the build route reads"""
    for compressed in (False, True):
        db = make_database(tmp_path / ("ann{gz}.fdb".format(gz=".gz" if compressed else "")))
        annotation_file = write(tmp_path / ("ann.tsv{gz}".format(gz=".gz" if compressed else "")),
                                ANNOTATION_ROWS, compressed)
        modify_module.ModifyTree(database=db, update_genomes=True).update_annotations(
            genomeid2taxid=annotation_file)
        rows = genome_rows(db)
        assert rows == {"GCF_000000001.1": (2, "refseq"), "GCA_000000002.1": (3, "genbank")}, \
            "compressed: {c}, rows: {r}".format(c=compressed, r=rows)


def test_every_row_of_a_compressed_modification_file_is_applied(tmp_path):
    """the file route of -mf adds the whole branch it describes, compressed or not: the family below the
    root and the species below the family, not only the last line of the file"""
    for compressed in (False, True):
        db = make_database(tmp_path / ("mod{gz}.fdb".format(gz=".gz" if compressed else "")))
        modfile = write(tmp_path / ("mod.tsv{gz}".format(gz=".gz" if compressed else "")),
                        MODIFICATION_ROWS, compressed)
        modify_object = modify_module.ModifyTree(database=db, mod_file=modfile, parent="root")
        modify_object.update_database()

        names = node_names(db)
        assert "Enterobacteriaceae" in names and "Escherichia coli" in names
        ids = ModifyFunctions(database=db).get_nodes(col=2)
        assert (ids["root"], ids["Enterobacteriaceae"]) in links(db), \
            "compressed: {c}, links: {l}".format(c=compressed, l=links(db))
        assert (ids["Enterobacteriaceae"], ids["Escherichia coli"]) in links(db)


def test_every_row_of_a_modification_file_without_levels_is_applied(tmp_path):
    """a modification file that names no level still describes a branch, and every row of it is applied
    (issue #142 kept the level empty, which had put the row handling out of the loop: only the last row
    of the file reached the database)"""
    db = make_database(tmp_path / "mod2col.fdb")
    modfile = write(tmp_path / "mod2col.tsv", [("parent", "child"),
                                               ("root", "Enterobacteriaceae"),
                                               ("Enterobacteriaceae", "Escherichia coli")])
    modify_module.ModifyTree(database=db, mod_file=modfile, parent="root").update_database()

    ids = ModifyFunctions(database=db).get_nodes(col=2)
    assert (ids["root"], ids["Enterobacteriaceae"]) in links(db), links(db)
    assert (ids["Enterobacteriaceae"], ids["Escherichia coli"]) in links(db), links(db)


def test_a_compressed_name_file_renames_nodes(tmp_path):
    """--update_names takes a file of old name, new name, and has to open it the same way"""
    rows = [("Shigella dysenteriae", "Shigella boydii")]
    for compressed in (False, True):
        db = make_database(tmp_path / ("names{gz}.fdb".format(gz=".gz" if compressed else "")))
        namefile = write(tmp_path / ("names.tsv{gz}".format(gz=".gz" if compressed else "")),
                         rows, compressed)
        modify_module.ModifyTree(database=db, update_node_names=True).update_node_names(namefile)
        names = node_names(db)
        assert "Shigella boydii" in names and "Shigella dysenteriae" not in names, \
            "compressed: {c}, names: {n}".format(c=compressed, n=names)


def test_a_modification_row_that_links_a_node_to_itself_is_skipped(tmp_path):
    """issue #162: the Brucella file opened with "Brucella melitensis	Brucella melitensis	species".
    A row whose parent and child are the same node is the one link shape no normal node may have, and
    it hides that node's subtree from the walk of the validation: the rest of the file is applied and
    the row is skipped, no self-link of a node other than the root is written"""
    db = make_database(tmp_path / "modself.fdb")
    modfile = write(tmp_path / "modself.tsv", [
        ("parent", "child", "level"),
        ("root", "Enterobacteriaceae", "family"),
        ("Enterobacteriaceae", "Enterobacteriaceae", "family"),
        ("Enterobacteriaceae", "Escherichia coli", "family"),
    ])
    modify_module.ModifyTree(database=db, mod_file=modfile, parent="root").update_database()
    conn = ModifyFunctions(database=db)
    assert conn.select("SELECT child FROM tree WHERE parent = child AND child != 1") == []
    ids = conn.get_nodes(col=2)
    assert (ids["root"], ids["Enterobacteriaceae"]) in links(db), links(db)
    assert (ids["Enterobacteriaceae"], ids["Escherichia coli"]) in links(db), links(db)
