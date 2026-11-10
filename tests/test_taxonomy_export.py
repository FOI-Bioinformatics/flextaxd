"""issue #134: the dmp export writes every link of the tree, also the levels that have no rank row.

The test taxonomy has the shape the pipeline really has: rank rows only for the official levels, and the
levels a custom taxonomy adds below species either have a rank row whose name is deliberately empty, have
no rank row at all, or have no rank_i. Reading the rank with an inner join dropped those links from
nodes.dmp, and the program reading the dump then lost that whole subtree without any message: kraken2
reported "Completed processing of 3 sequences" for four sequences and classified nothing below the node
that survived.
"""

import os
import random
import sys
import shutil
import subprocess

import pytest

from flextaxd.modules.WriteTaxonomy import WriteTaxonomy
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

NODES = [(1, "root"), (561, "Escherichia"), (562, "Escherichia coli"),
         (2001, "Escherichia coli subsp. A"), (2002, "Escherichia coli A strain-1"),
         (2003, "Escherichia coli A strain-2"), (3001, "Staphylococcus aureus"),
         (3002, "Staphylococcus aureus deep-1")]

"""rank_i 1 and 2 are the official levels, 3 is a rank row whose name is empty (the level below species
that the rank column is meant to leave empty), 4 is used by the tree and has no rank row at all"""
RANKS = [(1, "no rank"), (2, "species"), (3, "")]

"""(child, parent, rank_i): 2003 has a rank_i without a rank row, 3002 has no rank_i"""
TREE = [(1, 1, 1), (561, 1, 1), (562, 561, 2), (2001, 562, 3), (2002, 2001, 3),
        (2003, 2001, 4), (3001, 561, 2), (3002, 3001, None)]

"""the four leaves of the tree carry a sequence, the two levels below species included"""
LEAVES = {2001: "ACC2001", 2002: "ACC2002", 2003: "ACC2003", 3002: "ACC3002"}

SEED = 11
LENGTH = 700


def make_database(tmp_path):
    db = str(tmp_path / "custom.fdb")
    conn = ModifyFunctions(database=db)
    for taxid, name in NODES:
        conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(taxid, name))
    for rank_i, rank in RANKS:
        conn.query("INSERT INTO rank (rank_i,rank) VALUES (?,?)", insert_val=(rank_i, rank))
    for child, parent, rank_i in TREE:
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (?,?,?)", insert_val=(parent, child, rank_i))
    for taxid, accession in sorted(LEAVES.items()):
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, accession, "nt"))
    conn.commit()
    return db


def sequences():
    """a sequence per leaf, long enough for kraken2 (k=35, minimizer l=31) to index and match it"""
    random.seed(SEED)
    return {taxid: "".join(random.choice("ACGT") for _ in range(LENGTH)) for taxid in sorted(LEAVES)}


def dump(tmp_path, database, **kwargs):
    out = tmp_path / "dump"
    """the dump directory is created by the run that asks for the dump (the tmpdir of flextaxd-create)"""
    out.mkdir(exist_ok=True)
    writer = WriteTaxonomy(str(out), database=database, prefix="names,nodes", dbprogram="kraken2", **kwargs)
    writer.nodes()
    writer.names()
    return out


def node_lines(dump_dir):
    with open(str(dump_dir / "nodes.dmp")) as nodes:
        return [line.rstrip("\n").split("\t|\t") for line in nodes if line.strip()]


def test_every_link_of_the_tree_reaches_nodes_dmp(tmp_path):
    """the export writes one line per tree row, whatever rank the level has"""
    database = make_database(tmp_path)
    lines = node_lines(dump(tmp_path, database))

    assert len(lines) == len(TREE)
    assert sorted(int(line[0]) for line in lines) == sorted(child for child, _, _ in TREE)


def test_a_level_without_rank_keeps_an_empty_rank_field(tmp_path):
    """the ranks of the levels below species are empty on purpose: only the official levels have a rank row
    in the database, and the value that comes back for the others is NULL, which must not be printed as the
    text None (kraken2 numbers such levels itself, as S1, S2 ... of the species above them)"""
    database = make_database(tmp_path)
    lines = {int(line[0]): line for line in node_lines(dump(tmp_path, database))}

    for child, _, rank_i in TREE:
        assert child in lines
        rank_field = lines[child][2]
        assert rank_field != "None"
        if rank_i in (4, None):
            assert rank_field == ""
    """the levels that do have a rank keep it, and the empty rank row stays empty"""
    assert lines[562][2] == "species"
    assert lines[2001][2] == ""


def test_nodes_dmp_keeps_the_columns_kraken2_reads(tmp_path):
    """the trailing columns of the kraken2 dump are load bearing: the rank field has to stay in the third
    column, also for a line whose rank is empty (issue #136)"""
    database = make_database(tmp_path)
    lines = node_lines(dump(tmp_path, database))

    for line in lines:
        assert len(line) >= 5
        assert line[0].isdigit() and line[1].isdigit()


KRAKEN2 = shutil.which("kraken2-build")


@pytest.mark.skipif(not KRAKEN2, reason="kraken2 is not installed")
def test_kraken2_indexes_and_classifies_the_levels_without_rank(tmp_path):
    """the point of keeping the links: the sequences of the levels without an official rank are in the
    built database and reads find them, and kraken2 numbers those levels itself (S1, S2 below species)"""
    database = make_database(tmp_path)
    dump_dir = dump(tmp_path, database)
    db = tmp_path / "k2db"
    (db / "taxonomy").mkdir(parents=True)
    (db / "library").mkdir()
    shutil.copy(str(dump_dir / "nodes.dmp"), str(db / "taxonomy" / "nodes.dmp"))
    shutil.copy(str(dump_dir / "names.dmp"), str(db / "taxonomy" / "names.dmp"))

    seqs = sequences()
    with open(str(db / "library" / "library.fna"), "w") as library, \
            open(str(db / "library" / "prelim_map.txt"), "w") as prelim, \
            open(str(db / "reads.fa"), "w") as reads:
        for taxid, accession in sorted(LEAVES.items()):
            sequence = seqs[taxid]
            library.write(">{0}\n{1}\n".format(accession, "\n".join(
                sequence[start:start + 70] for start in range(0, LENGTH, 70))))
            """the taxid of every sequence travels in the map, as the library step of a build writes it"""
            prelim.write("TAXID\t{0}\t{1}\n".format(accession, taxid))
            reads.write(">read_{0}\n{1}\n".format(taxid, sequence[200:450]))

    build = subprocess.run([KRAKEN2, "--build", "--skip-maps", "--db", str(db), "--threads", "2"],
                           capture_output=True, text=True)
    assert "Completed processing of {0} sequences".format(len(LEAVES)) in build.stdout + build.stderr

    classify = subprocess.run([shutil.which("kraken2"), "--db", str(db), "--quick", str(db / "reads.fa")],
                              capture_output=True, text=True)
    hits = {line.split("\t")[1]: line.split("\t")[2]
            for line in classify.stdout.splitlines() if line.startswith("C")}
    for taxid in LEAVES:
        assert hits.get("read_{0}".format(taxid)) == str(taxid)

    inspect = subprocess.run([shutil.which("kraken2-inspect"), "--db", str(db), "--report-zero-counts",
                              "--threads", "2"], capture_output=True, text=True)
    report = inspect.stdout
    """the unnamed levels are kept with the names from names.dmp and the rank letters kraken2 gives them"""
    for taxid, name in NODES:
        if taxid in LEAVES:
            assert str(taxid) in report and name in report
    assert "\tS1\t" in report and "\tS2\t" in report

"""issue #136: the minimal layout and a program format exclude each other, and the run that asks for
both is stopped before a dump is written.

--dump_mini writes the taxonomy without the trailing columns of the dmp files, and those columns are
what kraken2 and krakenuniq read for the name of a taxid. Measured with kraken2 2.17.2 on the tree of
this file: the same four sequences get indexed either way, but the database built from a minimal dump
lists every taxid without a name and reads as an empty database. Such a dump is useless and was written
without a word about it.
"""


def test_a_minimal_dump_is_never_written_for_a_program(tmp_path):
    """the writer keeps both settings and refuses the combination instead of quietly dropping the
    columns the program reads"""
    from flextaxd.modules.WriteTaxonomy import DumpFormatError

    database = make_database(tmp_path)
    out = tmp_path / "mini"
    out.mkdir()
    writer = WriteTaxonomy(str(out), database=database, prefix="names,nodes", dbprogram="kraken2")
    with pytest.raises(DumpFormatError):
        writer.set_minimal()

    """the combination can also arrive through the command line (--dump_mini together with
    --dbprogram), and both writers refuse to write the file"""
    forced = WriteTaxonomy(str(out), database=database, prefix="names,nodes",
                           dbprogram="kraken2", minimal=True)
    with pytest.raises(DumpFormatError):
        forced.nodes()
    with pytest.raises(DumpFormatError):
        forced.names()
    assert not (out / "nodes.dmp").exists()

    """a program dump that is not minimal keeps the trailing columns the program reads, and the minimal
    dump of the same tree is shorter: the two layouts really do differ"""
    program_dir = tmp_path / "program"
    program_dir.mkdir()
    program = WriteTaxonomy(str(program_dir), database=database, prefix="names,nodes", dbprogram="kraken2")
    program.nodes()
    program.names()
    minimal_dir = tmp_path / "plain"
    minimal_dir.mkdir()
    plain = WriteTaxonomy(str(minimal_dir), database=database, prefix="names,nodes", minimal=True)
    plain.nodes()
    plain.names()

    wide = max(len(line.rstrip("\n").split("\t|\t"))
               for line in open(str(program_dir / "names.dmp")) if line.strip())
    narrow = max(len(line.rstrip("\n").split("\t"))
                 for line in open(str(minimal_dir / "names.dmp")) if line.strip())
    assert wide > narrow, "the minimal layout is not dropping the columns of the program format"
    assert len([line for line in open(str(minimal_dir / "nodes.dmp")) if line.strip()]) == len(TREE)


def test_the_dump_command_stops_before_writing_an_unusable_dump(tmp_path, monkeypatch):
    """the combination is refused where a clean message is possible: the command line for
    --dump_mini --dbprogram kraken2 ends without writing any dmp file"""
    from flextaxd import custom_taxonomy_databases

    database = make_database(tmp_path)
    out = tmp_path / "clidump"
    out.mkdir()
    monkeypatch.setattr(sys, "argv", ["flextaxd", "-db", database, "-o", str(out),
                                      "--dump_mini", "--dbprogram", "kraken2"])
    with pytest.raises(SystemExit):
        custom_taxonomy_databases.main()

    assert not (out / "nodes.dmp").exists()
    assert not (out / "names.dmp").exists()
