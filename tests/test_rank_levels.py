"""issue #138: a rank name gets one rank level, whatever import asks for it.

`DatabaseConnection.add_rank` inserted unconditionally, and only the importer object kept a name -> id
dictionary of its own, so every import run against the same database added its own copy of every rank it
named: gtdb.fdb (bug1/flextaxd_purge) holds `no rank` twice, rank_i 1 and 2, because two imports ran against
that file. The dump still prints the name of whichever row a link points at, but a level of the tree then
exists under several ids and anything reasoning about levels has to match names and can land on any of them.
The rank column has no unique constraint on purpose: the databases already on disk hold duplicates.
"""

import pytest

from flextaxd.modules.ReadTaxonomyCanSNPer import ReadTaxonomyCanSNPer
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

RANKS = [(1, "no rank"), (2, "genus"), (3, "species")]

TREE = """Enterococcus faecium;CC17;ST18
Enterococcus faecium;CC17;ST413
Enterococcus faecium;CC114;ST308
"""


def make_database(tmp_path, ranks=RANKS):
    db = str(tmp_path / "ranks.fdb")
    conn = ModifyFunctions(database=db)
    for rank_i, rank in ranks:
        conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(rank_i, rank))
    conn.commit()
    return db


def rank_rows(database):
    conn = ModifyFunctions(database=database)
    return conn.select("SELECT rank_i,rank FROM rank ORDER BY rank_i")


def rank_levels(database):
    """name -> id, and the number of rows says whether a name was added more than once"""
    rows = rank_rows(database)
    return {rank: rank_i for rank_i, rank in rows}, len(rows)


def make_tree(tmp_path):
    tree_file = tmp_path / "tree.txt"
    tree_file.write_text(TREE)
    return str(tree_file)


def test_an_import_uses_the_rank_level_it_already_has(tmp_path):
    """two imports of the same ranks must not leave two levels per name, they must point at the same row"""
    database = make_database(tmp_path)
    first = ModifyFunctions(database=database)
    asked = {name: first.add_rank(name) for name in ("no rank", "genus", "species", "strain")}
    first.commit()
    second = ModifyFunctions(database=database)
    added = {name: second.add_rank(name) for name in ("no rank", "genus", "species", "strain")}
    second.commit()

    levels, rows = rank_levels(database)
    assert rows == 4, "an import added its own copy of a rank level"
    assert sorted(levels) == ["genus", "no rank", "species", "strain"]
    """the ranks that were already in the database keep their ids, for both imports"""
    assert asked["no rank"] == added["no rank"] == 1
    assert asked["genus"] == added["genus"] == 2
    assert asked["species"] == added["species"] == 3
    assert asked["strain"] == added["strain"] and asked["strain"] > 3
    assert len(set(levels.values())) == len(levels), "two names share one rank level"


def test_an_import_into_a_database_with_ranks_adds_no_copy_of_them(tmp_path):
    """the CanSNPer import registers its ranks through the base importer, which used to insert a second
    `no rank` row next to the one the database already had - the duplicate of gtdb.fdb"""
    database = make_database(tmp_path)
    reader = ReadTaxonomyCanSNPer(taxonomy_file=make_tree(tmp_path), database=database, rank="species")
    reader.parse_taxonomy()

    levels, rows = rank_levels(database)
    assert rows == 3, "the import added copies of the rank levels it was given"
    assert sorted(levels) == ["genus", "no rank", "species"]
    assert levels["no rank"] == 1 and levels["species"] == 3

    conn = ModifyFunctions(database=database)
    links = conn.select("SELECT child,rank.rank FROM tree LEFT JOIN rank ON rank.rank_i = tree.rank_i")
    """the ranks that are named at all are the ones the database already had: the levels of the imported tree
    carry no rank (issue #142), and none of them landed on a rank of the backbone they were added to"""
    for _, rank in links:
        assert rank in (None, "no rank", "species")


def test_an_explicit_rank_level_is_written_to_the_level_column(tmp_path):
    """issue #139: a caller that names the level of a rank got an error instead of that level. The rank
    table has the column rank_i and no id, and the explicit id branch of add_rank wrote the value into a
    column called id, so the insert failed and the exception object came back as the level"""
    database = make_database(tmp_path, ranks=[])
    conn = ModifyFunctions(database=database)
    level = conn.add_rank("species", id=7)
    conn.commit()

    assert level == 7, "the caller did not get the level it asked for"
    levels, rows = rank_levels(database)
    assert rows == 1
    assert levels == {"species": 7}
    """the name keeps that level for later callers, whether they name it or not"""
    assert conn.add_rank("species") == 7
    assert conn.add_rank("genus") > 7


def test_a_level_without_a_rank_name_stays_one_level(tmp_path):
    """the empty rank name is the level the dump writes as an empty field (issue #134), and it is a name
    like any other: asking for it twice must not add a second one"""
    database = make_database(tmp_path, ranks=[])
    conn = ModifyFunctions(database=database)
    first = conn.add_rank("")
    second = conn.add_rank("")
    conn.commit()

    assert first == second
    levels, rows = rank_levels(database)
    assert rows == 1 and list(levels.values()) == [first]
    with pytest.raises(AttributeError):
        conn.add_rank(None)
