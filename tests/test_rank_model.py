"""the rank levels of a database, one model across the suite (issues #140, #141, #142).

A level of the tree is named by a rank or it is not named: NCBI writes "no rank" for the taxids it does not
place in a rank and that name is data worth keeping (238,697 such links in NCBI_base.db on disk), while a file
that carries no rank column at all - a CanSNPer tree, a modification line with two columns, a QIIME level code
outside the vocabulary - names nothing and the link of that level must carry no rank id. rank_i 1 used to be
the default of add_link, and rank_i 1 is the "no rank" level of every database derived from an NCBI taxdump,
so an unnamed level was silently named "no rank" instead of staying unnamed (#142), and the tree modification
code asked whether a rank was the fallback by comparing the rank *name* with the id 1, a test that never fired
(#140).
"""

import logging


from flextaxd.modules import ranks
from flextaxd.modules.ModifyTree import ModifyTree
from flextaxd.modules.ReadTaxonomyNCBI import ReadTaxonomyNCBI
from flextaxd.modules.ReadTaxonomyQIIME import ReadTaxonomyQIIME
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

"""the shape of an NCBI nodes.dmp (rank field followed by the empty gc mask field, as the real files have
it): NCBI's own "no rank" name, the official ranks, two rows whose rank field names nothing (an empty field,
and the text None as an export of this code once wrote it), and a rank name of a custom taxonomy"""
NODES_DMP = "".join(
	"{0}\t|\t{1}\t|\t{2}\t|\t\t|\n".format(child, parent, rank) for child, parent, rank in [
		(1, 1, "no rank"),
		(2, 1, "superkingdom"),
		(1224, 2, "no rank"),
		(561, 1224, "genus"),
		(562, 561, "species"),
		(14, 562, ""),
		(15, 562, "None"),
		(16, 562, "subsubspecies"),
	])

QIIME_FILE = (
	"GENOME1\td__Bacteria;g__Escherichia;s__Escherichia coli\n"
	"GENOME2\tg__Escherichia;z__Weird thing\n"
)


def make_database(tmp_path):
    return str(tmp_path / "ranks.fdb")


def tree_links(database):
    """(child, parent) -> rank name of the link, None when the link carries no rank level"""
    conn = ModifyFunctions(database=database)
    links = {}
    for child, parent, rank in conn.select("SELECT child,parent,rank FROM tree "
                                           "LEFT JOIN rank ON rank.rank_i = tree.rank_i"):
        links[(int(child), int(parent))] = rank
    return links


def rank_names(database):
    conn = ModifyFunctions(database=database)
    return {name for _, name in conn.select("SELECT rank_i,rank FROM rank")}


def test_what_ncbi_names_is_kept_and_what_ncbi_leaves_blank_stays_blank(tmp_path):
    database = make_database(tmp_path)
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(NODES_DMP)
    reader = ReadTaxonomyNCBI(taxonomy_file=str(nodes), database=database, force_multisource=False)
    reader.read_nodes(str(nodes))

    links = tree_links(database)
    assert links[(562, 561)] == "species"
    assert links[(1224, 2)] == "no rank", "the name NCBI gives its unnamed taxids is data, it is kept"
    assert links[(14, 562)] is None, "an empty rank field names nothing: the link carries no rank level"
    assert links[(15, 562)] is None, "the text None names nothing either"

    """a rank name outside the NCBI vocabulary is stored as it came in, and said out loud as custom"""
    assert "subsubspecies" in rank_names(database)
    assert links[(16, 562)] == "subsubspecies"
    summary = ranks.summary(ModifyFunctions(database=database))
    assert summary["custom"] == ["subsubspecies"]
    assert summary["unnamed"] == 2
    assert summary["dangling"] == 0, "no link may point at a rank level that is not in the rank table"




def test_the_ncbi_reader_can_be_built_without_the_multisource_option(tmp_path):
    """issue #144: the reader read force_multisource out of the keyword arguments it required, and only
    the command line of flextaxd-create passes that option, so building the reader any other way raised
    TypeError before a single line was read"""
    database = make_database(tmp_path)
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(NODES_DMP)
    reader = ReadTaxonomyNCBI(taxonomy_file=str(nodes), database=database)
    assert reader.force_multisource is False, "the option is not required and defaults to off"

    """and it reads the same way it always did"""
    reader.read_nodes(str(nodes))
    assert tree_links(database)[(562, 561)] == "species"
    assert tree_links(database)[(1224, 2)] == "no rank"

    """a caller that does pass the option keeps it"""
    kept = ReadTaxonomyNCBI(taxonomy_file=str(nodes), database=database, force_multisource=True)
    assert kept.force_multisource is True

def test_the_vocabulary_covers_the_ranks_of_the_databases_on_disk():
    """the 44 rank names of the NCBI taxdump in NCBI_base.db, plus the domain level the QIIME/GTDB reader
    uses: none of them may be reported as custom, and the ranks a user adds themselves may"""
    assert ranks.custom_ranks(ranks.OBSERVED) == []
    assert ranks.is_defined("domain") and ranks.is_defined("subspecies") and ranks.is_defined("forma specialis")
    assert not ranks.is_defined("no rank"), "the NCBI name for an unranked taxid is not a rank of the tree"
    assert ranks.custom_ranks(["subsubspecies", "clade", "no rank", ""]) == ["subsubspecies"]


def test_a_link_takes_no_rank_id_by_default(tmp_path):
    """rank_i 1 was the default of add_link and is the "no rank" level of an NCBI derived database (#142)"""
    database = make_database(tmp_path)
    conn = ModifyFunctions(database=database)
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(1, "no rank"))
    conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(2, "species"))
    parent = conn.add_node("Escherichia")
    child = conn.add_node("Escherichia coli")
    conn.add_link(child=child, parent=parent)
    conn.commit()

    links = tree_links(database)
    assert links[(child, parent)] is None, "an unnamed link lands on a rank level of its own account"


def test_a_qiime_level_outside_the_vocabulary_is_reported_and_not_fatal(tmp_path, caplog):
    """a level code the QIIME/GTDB format does not define used to reach add_rank as None, which raised
    "Rank cannot be None" in the middle of an import (#141)"""
    database = make_database(tmp_path)
    tree_file = tmp_path / "taxonomy.txt"
    tree_file.write_text(QIIME_FILE)

    with caplog.at_level(logging.WARNING, logger="flextaxd.modules.ReadTaxonomy"):
        reader = ReadTaxonomyQIIME(taxonomy_file=str(tree_file), database=database)
        reader.set_qiime(True)
        reader.parse_taxonomy()

    assert "not a QIIME level code" in caplog.text
    conn = ModifyFunctions(database=database)
    nodes = {name: taxid for taxid, name in conn.select("SELECT id,name FROM nodes")}
    assert "Weird thing" in nodes, "the node is still added"
    unnamed_children = {child for child, in conn.select("SELECT child FROM tree WHERE rank_i IS NULL")}
    assert nodes["Weird thing"] in unnamed_children, "the level it could not name stays without a rank"
    """the levels the format does define keep their names"""
    assert "species" in rank_names(database) and "genus" in rank_names(database)


def test_no_rank_is_not_a_level_of_the_tree(tmp_path):
    """keep_levels decides which level a re-added node keeps from the node above it: "no rank" and a link
    with no rank at all are not levels (#140)"""
    database = make_database(tmp_path)
    conn = ModifyFunctions(database=database)
    for rank_i, name in ((1, "no rank"), (2, "genus"), (3, "species")):
        conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(rank_i, name))
    conn.commit()

    modifier = ModifyTree(database=database, clean_database=True)
    links = [(11, 1, 1), (12, 2, None), (13, 3, 3), (14, 4, 2)]
    levels = modifier.keep_levels(links)

    assert levels == {3: 3, 4: 2}, "only the parents that really carry a level are kept"
    assert ranks.is_no_rank(None) and ranks.is_no_rank("no rank") and ranks.is_no_rank("")
    assert not ranks.is_no_rank("genus")
