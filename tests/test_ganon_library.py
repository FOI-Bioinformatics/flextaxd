"""issue #133 and #146: a build with --dbprogram ganon died before any work was done.

create_databases.py hands one argument list to whichever program was selected, and CreateGanonDB
accepted neither tmpdir, create_lib, skip_genomes nor multifiles, so every --dbprogram ganon run ended
with "TypeError: unexpected keyword argument". The library step of ganon also never saw the packed
files and the build listed its library files by hand over the process numbers, which left the packed
sequences out of what ganon reads. ganon itself is not installed and not called here: the assertions
are about the arguments the class takes, the library and map it writes, and the start method the parallel
header pass asks for.
"""

import glob
import multiprocessing
import os
import shutil
import sys

import pytest

from subprocess import CalledProcessError

from flextaxd.modules.CreateGanonDB import CreateGanonDB
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

TESTDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testdata")
GTDB_GENOMES = os.path.join(TESTDATA, "gtdb", "genomes")
NT_PACKED = os.path.join(TESTDATA, "nt", "nt.fna.gz")

"""the four genomes of the shipped folder that have an annotation in the tree below"""
GENOME_TAXIDS = {
    "GCF_000005845.2": 562,
    "GCF_000013465.1": 1280,
    "GCA_000016685.1": 1280,
    "GCF_000009065.1": 1281,
}

"""the sequence level annotation of the packed file, as an accession2taxid file gives it"""
NT_TAXIDS = {
    "NC_000913.3": 511145,
    "NC_007795.1": 367579,
    "NC_004470.1": 1841072,
    "CP000250.1": 562,
}

UNANNOTATED_SEQUENCE = "LC_test_1"


def make_database(tmp_path):
    """A small database holding the genome annotations of the shipped folder and the sequence level
    annotations of the packed file"""
    db = str(tmp_path / "ganon.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    for taxid in sorted(set(list(GENOME_TAXIDS.values()) + list(NT_TAXIDS.values()))):
        conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(taxid, "taxid" + str(taxid)))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,?,1)", insert_val=(taxid,))
    for genome, taxid in sorted(GENOME_TAXIDS.items()):
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, genome, "refseq"))
    for seqid, taxid in sorted(NT_TAXIDS.items()):
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, seqid, "nt"))
    conn.commit()
    return db


def genome_folder(tmp_path):
    """The genome files of the build, as create_databases hands them over: a dictionary of genome id
    to path. The id of an official file is its accession, the `_genomic` suffix of the file name is not
    part of it (ProcessDirectory.is_gcf_gca)"""
    names = {}
    for path in sorted(os.listdir(GTDB_GENOMES)):
        if not path.endswith(".fna"):
            continue
        name = path[:-len(".fna")]
        accession = name.split("_genomic")[0]
        key = accession if accession.startswith(("GCF_", "GCA_")) else name
        names[key] = os.path.join(GTDB_GENOMES, path)
    return names


def map_lines(path):
    lines = {}
    if not os.path.exists(str(path)):
        return lines
    with open(str(path)) as mapfile:
        for row in mapfile:
            columns = row.rstrip("\n").split("\t")
            if len(columns) >= 3:
                lines[columns[0]] = columns[2]
    return lines


def ganon_build(tmp_path, database, genome_names, **kwargs):
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    ganondb = tmp_path / "ganondb"
    kwargs.setdefault("processes", 1)
    return CreateGanonDB(database, str(ganondb), genome_names, str(out), **kwargs), ganondb, out


@pytest.fixture
def packed(tmp_path):
    """the packed file of the test, copied so the shipped testdata stays untouched"""
    target = tmp_path / "nt.fna.gz"
    shutil.copy(NT_PACKED, str(target))
    return str(target)


def test_ganon_takes_the_arguments_of_the_shared_build_interface(tmp_path, packed):
    """The argument list create_databases.py builds for the selected program has to be accepted by
    ganon too: tmpdir, create_lib, skip_genomes and multifiles were refused outright (#133)."""
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(
        tmp_path, database, genome_folder(tmp_path),
        create_db=True, limit=0, dbprogram="ganon", params="", skip=False,
        processes=1, build_processes=2, debug=False, verbose=False,
        tmpdir=str(tmp_path / "tmp"), create_lib=False,
        skip_genomes=["GCF_000013465.1"], multifiles=[packed])

    assert build.tmpdir == str(tmp_path / "tmp")
    assert build.skip_genomes == ["GCF_000013465.1"]
    assert build.packed_files == [packed]
    assert build.build_processes == 2


def test_packed_sequences_reach_the_ganon_library_and_its_map(tmp_path, packed):
    """The packed file is not split into genome files: its sequences go into the library as they stand
    and their taxid (plus the length ganon asks for) goes into the seq-info map."""
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(tmp_path, database, genome_folder(tmp_path))
    build.create_library_from_files(multifiles=[packed])

    """the packed sequences are in the library with the header the archive has, no ganon suffix and no
    kraken tag, and unannotated sequences are not written at all"""
    packed_library = ganondb / ".tmppacked.fasta"
    assert packed_library.exists()
    headers = [line.rstrip("\n") for line in open(str(packed_library)) if line.startswith(">")]
    assert sorted(headers) == [">" + seqid for seqid in sorted(NT_TAXIDS)]
    assert UNANNOTATED_SEQUENCE not in "".join(headers)

    """the map line of a packed sequence is its own id, its length and its taxid; the genome sequences
    of the same run keep the unique id ganon wants (sequence id plus genome id)"""
    taxids = map_lines(build.seqid2taxid)
    for seqid, taxid in NT_TAXIDS.items():
        assert taxids[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in taxids
    assert any(header.endswith("_" + genome) for header in taxids for genome in GENOME_TAXIDS)

    """ganon reads the length of every sequence out of that map, which the kraken style map does not
    carry: the line of a packed sequence names the length of the sequence as it stands in the library"""
    lengths = {}
    for row in open(str(build.seqid2taxid)):
        columns = row.rstrip("\n").split("\t")
        if len(columns) >= 2 and columns[0] in NT_TAXIDS:
            lengths[columns[0]] = int(columns[1])
    assert sorted(lengths) == sorted(NT_TAXIDS)
    assert all(length > 0 for length in lengths.values()), "the map carries no sequence length"

    """the build hands ganon every library file of the run: the per process files and the packed file,
    whichever name they ended up with (#133)"""
    library_files = [os.path.basename(path) for path in glob.glob(str(ganondb / ".tmp*.fasta*"))]
    assert ".tmppacked.fasta" in library_files
    assert any(name != ".tmppacked.fasta" for name in library_files), \
        "the library files of the genome processes are not part of what the build reads: {f}".format(f=library_files)


def test_genomes_on_the_skip_list_are_left_out_of_the_ganon_library(tmp_path, packed):
    """--skip_genomes meant nothing to the ganon class: skipped genomes still entered the library."""
    genomes = genome_folder(tmp_path)
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(tmp_path, database, genomes,
                                      skip_genomes=["GCF_000013465.1"])
    build.create_library_from_files(multifiles=[])

    taxids = map_lines(build.seqid2taxid)
    assert not any(header.endswith("_GCF_000013465.1") for header in taxids)
    assert any(header.endswith("_" + genome) for genome in genomes
               if genome != "GCF_000013465.1" for header in taxids)


def test_the_ganon_header_pass_names_the_start_method_it_needs(tmp_path, monkeypatch, packed):
    """The children of the header pass inherit the database connection and the genome paths of the build
    instead of having them pickled over, so the pass cannot rely on the default start method: python 3.14
    made forkserver the default, and that pickles the target and dies on the sqlite connection (issue #146,
    the failure #130 fixed for kraken). The pass asks for fork itself, and the library shows the children
    really ran."""
    requested = []
    real_get_context = multiprocessing.get_context

    def spy(name=None):
        requested.append(name)
        return real_get_context(name)

    monkeypatch.setattr(multiprocessing, "get_context", spy)
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(tmp_path, database, genome_folder(tmp_path), processes=2)
    build.create_library_from_files(multifiles=[packed])

    assert "fork" in requested, "the ganon header pass left the start method to the interpreter"
    taxids = map_lines(build.seqid2taxid)
    for genome in GENOME_TAXIDS:
        assert any(header.endswith("_" + genome) for header in taxids), \
            "{0} did not get a map line from the header pass".format(genome)
    library_files = sorted(os.path.basename(path) for path in glob.glob(str(ganondb / ".tmp*.fasta")))
    assert ".tmppacked.fasta" in library_files


"""the two files of the taxdump that create_database copies into the database directory"""
def write_taxonomy(outdir):
    outdir.mkdir(exist_ok=True)
    (outdir / "1_nodes.dmp").write_text("1\t|\troot\t|\tno rank\t|\n")
    (outdir / "1_names.dmp").write_text("1\t|\troot\t|\tscientific name\t|\n")


def fake_ganon(tmp_path, monkeypatch, behaviour="ok"):
    """a ganon on the PATH standing in for the real one, which is not installed here: it writes the marker
    file where the -d argument points, or fails the way the issue measured (exit 3 with a message on
    stderr). A build that stops is what the old os.system call threw away."""
    bin = tmp_path / "bin"
    bin.mkdir(exist_ok=True)
    script = bin / "ganon"
    script.write_text("""#!{python}
import os, sys
args = sys.argv[1:]
if not args or args[0] != "build":
    sys.exit(0)
if "{behaviour}" == "fail":
    print("ganon: the genome files could not be read", file=sys.stderr)
    sys.exit(3)
db = args[args.index("-d") + 1]
open(os.path.join(db, "fake.ganondb"), "w").write("built")
""".format(python=sys.executable, behaviour=behaviour))
    os.chmod(str(script), 0o755)
    monkeypatch.setenv("PATH", str(bin) + os.pathsep + os.environ.get("PATH", ""))


def ganon_with_library(tmp_path, packed):
    """a build that has written its library and its taxonomy dump, as create_databases has it when it
    reaches the build step"""
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(tmp_path, database, genome_folder(tmp_path), multifiles=[packed])
    build.create_library_from_files(multifiles=[packed])
    write_taxonomy(out)
    return build, ganondb, out


def test_a_ganon_that_is_not_installed_stops_the_run(tmp_path, packed, monkeypatch):
    """issue #148: the build ran through os.system and its exit status was thrown away, so on a machine
    without ganon the run printed "sh: 1: ganon: not found", deleted every library file it had written and
    still ended with the database reported as built. The CalledProcessError is what create_databases
    already turns into "no ganon database was created" and a failing exit status (issue #132)."""
    build, ganondb, out = ganon_with_library(tmp_path, packed)
    (tmp_path / "emptybin").mkdir(exist_ok=True)
    monkeypatch.setenv("PATH", str(tmp_path / "emptybin"))
    with pytest.raises(CalledProcessError) as failure:
        build.create_database(str(out))
    assert failure.value.returncode == 127
    assert glob.glob(str(ganondb / ".tmp*.fasta*")), "the library of the run was deleted on the way out"


def test_a_failing_ganon_build_stops_the_run(tmp_path, packed, monkeypatch):
    """a ganon that stops with a message is the same failure: the run has no database and says so, with
    what ganon wrote on stderr"""
    build, ganondb, out = ganon_with_library(tmp_path, packed)
    fake_ganon(tmp_path, monkeypatch, behaviour="fail")
    with pytest.raises(CalledProcessError) as failure:
        build.create_database(str(out))
    assert failure.value.returncode == 3
    assert glob.glob(str(ganondb / ".tmp*.fasta*")), "the library of the run was deleted on the way out"


def test_a_ganon_build_keeps_the_library_it_wrote(tmp_path, packed, monkeypatch):
    """the library files are the library of this build: the next database should be able to use them
    (T6), so they stay where they are and the taxonomy reaches the database directory"""
    build, ganondb, out = ganon_with_library(tmp_path, packed)
    fake_ganon(tmp_path, monkeypatch)
    build.create_database(str(out))

    assert (ganondb / "fake.ganondb").exists(), "the build was not run over the library"
    library_files = sorted(os.path.basename(path) for path in glob.glob(str(ganondb / ".tmp*.fasta*")))
    assert ".tmppacked.fasta" in library_files, "the library of the run is gone: {0}".format(library_files)
    taxonomy = ganondb / "taxonomy"
    assert (taxonomy / "nodes.dmp").exists() and (taxonomy / "names.dmp").exists()


def test_a_ganon_build_without_any_library_is_not_started(tmp_path, monkeypatch):
    """the old code logged the empty -i argument and asked ganon to build an empty database anyway"""
    database = make_database(tmp_path)
    build, ganondb, out = ganon_build(tmp_path, database, genome_folder(tmp_path))
    ganondb.mkdir(exist_ok=True)
    write_taxonomy(out)
    fake_ganon(tmp_path, monkeypatch)
    with pytest.raises(CalledProcessError):
        build.create_database(str(out))
    assert not (ganondb / "fake.ganondb").exists(), "ganon was started over an empty library"
