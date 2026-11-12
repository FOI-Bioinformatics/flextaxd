"""T3: a build given a genome folder and an already packed file puts both into library.fna and gives
every sequence that has an annotation its line in the taxid map, without touching the packed file.

The build is run through the command line of create_databases.py, because that is where the packed
file used to be written into the library by itself (setting skip=True, issue #114) and where the taxid
map used to get lines for genome files only (issue #115). kraken2 itself is not called: the assertions
are about the library and the map the build leaves behind.
"""

import builtins
import glob
import gzip
import hashlib
import multiprocessing
import os
import random
import re
import shutil
import sqlite3
import subprocess
import sys

import pytest

from flextaxd import create_databases
from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

"""kraken2 has to be installed for the two tests that ask the program what it made of the library
(issue #135); the other tests of this file only look at the files the build leaves behind"""
KRAKEN2 = shutil.which("kraken2-build")

TESTDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testdata")
GTDB_GENOMES = os.path.join(TESTDATA, "gtdb", "genomes")
NT_PACKED = os.path.join(TESTDATA, "nt", "nt.fna.gz")

"""The six genome annotations of the GTDB test taxonomy (tests/testdata/gtdb/bac120_taxonomy.tsv), on
the NCBI taxids their species names give; the last two have no sequence file on disk"""
GENOME_TAXIDS = {
    "GCF_000005845.2": 562,
    "GCF_000013465.1": 1280,
    "GCA_000016685.1": 1280,
    "GCF_000009065.1": 1281,
    "GCF_000001415.1": 1423,
    "GCF_000009605.1": 28901,
}

"""The sequence level annotation of the packed file, as an accession2taxid file gives it"""
NT_TAXIDS = {
    "NC_000913.3": 511145,
    "NC_007795.1": 367579,
    "NC_004470.1": 1841072,
    "CP000250.1": 562,
}

UNANNOTATED_SEQUENCE = "LC_test_1"

NODE_NAMES = [(562, "Escherichia coli"),
              (1280, "Staphylococcus aureus"),
              (1281, "Staphylococcus epidermidis"),
              (1423, "Bacillus subtilis"),
              (28901, "Salmonella enterica"),
              (511145, "Escherichia coli O157:H7 str. TW13578"),
              (367579, "Staphylococcus aureus subsp. aureus NCTC 8325"),
              (1841072, "Staphylococcus epidermidis ATCC 12228")]


def make_database(tmp_path):
    """A database holding the six GTDB genome annotations and the four sequence level annotations of
    the packed file, every one of them on a node that exists in the tree"""
    db = str(tmp_path / "build.fdb")
    conn = ModifyFunctions(database=db)
    conn.query("INSERT INTO nodes (id,name) VALUES (1,'root')")
    conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,1,1)")
    for taxid, name in NODE_NAMES:
        conn.query("INSERT INTO nodes (id,name) VALUES (?,?)", insert_val=(taxid, name))
        conn.query("INSERT INTO tree (parent,child,rank_i) VALUES (1,?,1)", insert_val=(taxid,))
    for genome, taxid in sorted(GENOME_TAXIDS.items()):
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, genome, "refseq"))
    for seqid, taxid in sorted(NT_TAXIDS.items()):
        conn.query("INSERT INTO genomes (id,genome,reference) VALUES (?,?,?)",
                   insert_val=(taxid, seqid, "nt"))
    conn.commit()
    return db


def copy_inputs(tmp_path):
    """Copy the shipped test data into the test directory: the build must not touch the packed file and
    the shipped testdata stays as it is for the other tests"""
    genomes = tmp_path / "genomes"
    shutil.copytree(GTDB_GENOMES, str(genomes))
    packed = tmp_path / "nt.fna.gz"
    shutil.copy(NT_PACKED, str(packed))
    return genomes, packed


def run_build(tmp_path, monkeypatch, database, genomes_path=False, nt_source=False,
              multifile_prefix=False, create_db=True, answers=("n",)):
    """Run one flextaxd-create build and return the directories it wrote. The questions a build asks
    here are answered from `answers` (default: "n" to the question whether to download the two genomes
    that have no file on disk), everything else runs without interaction. --keep keeps the temp
    directory, which holds the dump of the map the build also writes."""
    work = tmp_path / "build"
    work.mkdir(exist_ok=True)
    out, tmpdir, dbname = work / "out", work / "tmp", work / "k2db"
    argv = ["flextaxd-create", "-db", database, "--db_name", str(dbname),
            "-o", str(out), "--tmpdir", str(tmpdir), "--keep"]
    if genomes_path:
        argv += ["--genomes_path", str(genomes_path)]
    if nt_source:
        argv += ["-nt", str(nt_source)]
    if multifile_prefix:
        argv += ["--multifile_prefix", multifile_prefix]
    argv.append("--create_db" if create_db else "--create_lib")
    answers = list(answers)

    def answer(prompt=""):
        return answers.pop(0) if answers else "n"

    monkeypatch.chdir(work)
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(builtins, "input", answer)
    create_databases.main()
    return out, tmpdir, dbname


def library_headers(dbname):
    with open(str(dbname / "library" / "library.fna")) as library:
        return [line.rstrip("\n") for line in library if line.startswith(">")]


def archive_headers(packed):
    with gzip.open(str(packed), "rt") as archive:
        return [line.rstrip("\n") for line in archive if line.startswith(">")]


def sequence_id(header):
    '''The accession of a library sequence: the first word of its header, without the kraken tag the
        genome files carry'''
    return header[1:].split()[0].split("|")[0]


def map_lines(*paths):
    """The taxid map of the build as a dictionary: sequence or genome id to taxid. Every map file the
    build writes is read (the dump of the database names the genomes, the map of the build names the
    sequences as they stand in library.fna), a later file wins over an earlier one for the same id."""
    lines = {}
    for path in paths:
        if not os.path.exists(str(path)):
            continue
        with open(str(path)) as mapfile:
            for row in mapfile:
                columns = row.rstrip("\n").split("\t")
                if len(columns) < 2:
                    continue
                key = columns[-2] if columns[0] == "TAXID" else columns[0]
                lines[key] = columns[-1]
    return lines


def build_map(out, tmpdir):
    return map_lines(out / "seqid2taxid.map", tmpdir / "seqid2taxid.map")


def md5(path):
    with open(str(path), "rb") as handle:
        return hashlib.md5(handle.read()).hexdigest()


def annotation_nodes_missing(database):
    conn = sqlite3.connect(str(database))
    return set(row[0] for row in conn.execute(
        "SELECT DISTINCT id FROM genomes WHERE id NOT IN (SELECT id FROM nodes)"))


def test_build_adds_genome_folder_and_packed_file_to_the_library(tmp_path, monkeypatch):
    """The Done when: case of T3. The library holds the four GTDB genome files and the five sequences
    of nt/nt.fna.gz, the map holds the six genome ids of the database plus the four annotated sequences
    of the packed file, and LC_test_1, which has no annotation, is only in the report."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    packed_before = md5(packed)

    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed)

    """the four genome files of the folder (two of the six annotated genomes have no file) plus the
    five sequences of the packed file; unannotated_genome.fna has no annotation and stays out"""
    headers = library_headers(dbname)
    assert sorted(sequence_id(header) for header in headers) == sorted(
        ["BA000036.1", "NC_000913.3", "NC_004470.1", "NC_007795.1"] +
        ["NC_000913.3", "NC_007795.1", "NC_004470.1", "CP000250.1", UNANNOTATED_SEQUENCE])

    """the genome files get the kraken header with their taxid, the packed sequences keep the header
    the archive has: their taxid travels in the map instead (the packed file is never split up)"""
    tagged = [header for header in headers if "kraken:taxid" in header]
    assert len(tagged) == 4
    """(the header line of a genome file ends with two spaces the build leaves there, hence the strip)"""
    assert sorted(header.split("|")[2].strip() for header in tagged) == sorted(
        str(GENOME_TAXIDS[genome]) for genome in
        ("GCF_000005845.2", "GCF_000013465.1", "GCA_000016685.1", "GCF_000009065.1"))
    assert [header for header in headers if "kraken:taxid" not in header] == archive_headers(packed)

    taxids = build_map(out, tmpdir)
    for genome, taxid in GENOME_TAXIDS.items():
        assert taxids[genome] == str(taxid)
    for seqid, taxid in NT_TAXIDS.items():
        assert taxids[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in taxids
    assert len(taxids) == len(GENOME_TAXIDS) + len(NT_TAXIDS)

    """the library step writes its own map, keyed by the ids exactly as library.fna holds them: the four
    sequences the archive brought in, and no sequence that has no annotation"""
    assert map_lines(out / "seqid2taxid.map") == {seqid: str(taxid) for seqid, taxid in NT_TAXIDS.items()}

    """the packed file was read, not touched"""
    assert md5(packed) == packed_before

    """the sequence without annotation is counted in the report and listed next to the input"""
    stored = ModifyFunctions(database=database)
    listed = str(packed) + ".unannotated"
    assert open(listed).read().split() == [UNANNOTATED_SEQUENCE]
    assert stored.get_report("unannotated_sequence", str(packed)) == [
        (str(packed), "total\t1\t" + listed)]
    report = "\n".join(stored.report_lines())
    assert "Library sequences without taxid: 1" in report
    assert "Sequence files without genome annotation: 1" in report
    """the disk report counts the two GTDB genome annotations that have no file of their own; the four
    sequence level annotations are not missing genomes, their sequences came in with the packed file and
    their taxids travel in the map, so they are reported on their own line (issue #131)"""
    assert stored.get_report("missing_on_disk") == [(str(genomes), "total\t2")]
    assert "Genome annotations without sequence file: 2" in report
    assert stored.get_report("sequence_annotation") == [(str(genomes), "total\t4")]
    assert stored.get_meta("sequence_annotations_on_disk") == "4"
    assert "Sequence level annotations (packed file, no genome file expected): 4" in report
    assert annotation_nodes_missing(database) == set()


def test_map_of_the_packed_file_reaches_the_kraken_library(tmp_path, monkeypatch):
    """kraken2 reads the taxids of a library from library/prelim_map.txt beside it: the sequences of
    the packed file are in the library with their own headers, so their lines have to be in that map."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed)

    prelim = map_lines(dbname / "library" / "prelim_map.txt")
    for seqid, taxid in NT_TAXIDS.items():
        assert prelim[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in prelim


@pytest.mark.skipif(not KRAKEN2, reason="kraken2 is not installed")
def test_the_program_of_the_build_indexes_the_library(tmp_path, monkeypatch, capfd):
    """issue #135: a build that writes a library and then produces nothing used to pass. What kraken2
    reports for the library is what the built database can find, so the build test reads it back."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed)

    captured = capfd.readouterr()
    indexed = re.search(r"Completed processing of (\d+) sequences", captured.out + captured.err)
    assert indexed, "the kraken2 step of the build reported nothing about the library it was given"

    """every sequence of the library that has a taxid has to be indexed: a genome file sequence carries
    its taxid in the header, a packed sequence finds it in prelim_map.txt beside the library. The one
    sequence of the packed file without annotation is the only one the program cannot place"""
    headers = library_headers(dbname)
    prelim = map_lines(dbname / "library" / "prelim_map.txt")
    annotated = [header for header in headers
                 if "kraken:taxid" in header or sequence_id(header) in prelim]
    assert len(annotated) == len(headers) - 1, "only the unannotated sequence is left out of the taxids"
    assert int(indexed.group(1)) == len(annotated), "the program indexed {} of the {} placed sequences".format(
        indexed.group(1), len(annotated))
    assert (dbname / "hash.k2d").exists(), "a build that indexed nothing leaves no database"


@pytest.mark.skipif(not KRAKEN2, reason="kraken2 is not installed")
def test_reads_of_a_library_of_real_sequence_find_their_genome(tmp_path, monkeypatch):
    """issue #135: the shipped test genomes are 160 bp of the repeat ACGTACGT, which kraken2 (k=35,
    minimizer l=31) indexes but never matches, so a read check cannot be made on them. Two genomes get
    a random sequence of real length here, and their reads have to come back on their taxids."""
    random.seed(7)
    sequences = {genome: "".join(random.choice("ACGT") for _ in range(700))
                 for genome in ("GCF_000005845.2", "GCF_000013465.1")}
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    for genome, sequence in sequences.items():
        (genomes / (genome + ".fna")).write_text(">{0}\n{1}\n".format(
            genome, "\n".join(sequence[start:start + 70] for start in range(0, 700, 70))))

    database = make_database(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database, genomes_path=genomes)

    reads = tmp_path / "reads.fa"
    reads.write_text("".join(">read_{0}\n{1}\n".format(genome, sequence[200:450])
                             for genome, sequence in sorted(sequences.items())))
    classify = subprocess.run([shutil.which("kraken2"), "--db", str(dbname), "--quick", str(reads)],
                              capture_output=True, text=True)
    hits = {line.split("\t")[1]: line.split("\t")[2]
            for line in classify.stdout.splitlines() if line.startswith("C")}
    for genome in sequences:
        assert hits.get("read_{0}".format(genome)) == str(GENOME_TAXIDS[genome]), \
            "the read of {0} did not come back on taxid {1}".format(genome, GENOME_TAXIDS[genome])


def test_library_only_build_also_maps_the_packed_file(tmp_path, monkeypatch):
    """--create_lib writes the library without calling the program: its map then also names the
    sequences of the genome files, while the packed file keeps the map lines of its accessions."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed, create_db=False)

    assert sorted(sequence_id(header) for header in library_headers(dbname)) == sorted(
        ["BA000036.1", "NC_000913.3", "NC_004470.1", "NC_007795.1"] +
        ["NC_000913.3", "NC_007795.1", "NC_004470.1", "CP000250.1", UNANNOTATED_SEQUENCE])
    taxids = build_map(out, tmpdir)
    for seqid, taxid in NT_TAXIDS.items():
        assert taxids[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in taxids
    """the genome sequences are in the map of the build too, with the taxid of their genome"""
    assert map_lines(out / "seqid2taxid.map")["BA000036.1"] == str(GENOME_TAXIDS["GCA_000016685.1"])


def test_packed_file_without_a_genome_folder_still_gets_its_map(tmp_path, monkeypatch):
    """A build with only the packed file keeps working (it used to cat the archive into library.fna and
    write nothing to the map), the genome folder is simply not used."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database, nt_source=packed)

    assert [sequence_id(header) for header in library_headers(dbname)] == [
        "NC_000913.3", "NC_007795.1", "NC_004470.1", "CP000250.1", UNANNOTATED_SEQUENCE]
    taxids = build_map(out, tmpdir)
    for seqid, taxid in NT_TAXIDS.items():
        assert taxids[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in taxids
    assert open(str(packed) + ".unannotated").read().split() == [UNANNOTATED_SEQUENCE]


def test_a_library_losing_many_of_its_sequences_says_why(tmp_path, monkeypatch, caplog):
    """issue #71: when the annotation covers only part of the packed library, the build has to say what
    that usually is - an annotation file that describes only nucl_gb accessions - instead of leaving the
    user to count the .unannotated list and guess"""
    import logging

    database = make_database(tmp_path)
    packed = tmp_path / "many.fna"
    with open(str(packed), "w") as library_file:
        for accession in NT_TAXIDS:
            print(">{accession} annotated sequence".format(accession=accession), file=library_file)
            print("ACGT" * 30, file=library_file)
        for number in range(90):
            print(">CPX{number:04d}.1 not in the annotation file".format(number=number), file=library_file)
            print("ACGT" * 30, file=library_file)

    with caplog.at_level(logging.WARNING):
        run_build(tmp_path, monkeypatch, database, nt_source=packed)

    warnings = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    advice = [warning for warning in warnings if "nucl_wgs.accession2taxid" in warning]
    assert advice, warnings
    assert "90 of 94" in advice[0], advice[0]
    assert len(open(str(packed) + ".unannotated").read().split()) == 90


def test_a_library_with_a_stray_unannotated_sequence_gives_no_advice(tmp_path, monkeypatch, caplog):
    """the shipped fixture has one unannotated sequence in five, which is a sequence, not an annotation
    file that does not cover the library: the advice would be noise here"""
    import logging

    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    with caplog.at_level(logging.WARNING):
        run_build(tmp_path, monkeypatch, database, nt_source=packed)

    warnings = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    assert not [warning for warning in warnings if "nucl_wgs.accession2taxid" in warning], warnings


def test_merge_answer_adds_the_packed_file_to_an_existing_library(tmp_path, monkeypatch):
    """The interactive merge answer keeps the sequences already in the library and adds the packed file
    to it; its taxid map lines used to be missing entirely, the file was only cat'd onto the library."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database, genomes_path=genomes)
    first = library_headers(dbname)
    assert len(first) == 4
    assert not os.path.exists(str(packed) + ".unannotated")

    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database, genomes_path=genomes,
                                    nt_source=packed, answers=["n", "m"])
    headers = library_headers(dbname)
    assert headers[:4] == first
    assert [sequence_id(header) for header in headers[4:]] == [
        "NC_000913.3", "NC_007795.1", "NC_004470.1", "CP000250.1", UNANNOTATED_SEQUENCE]
    taxids = map_lines(out / "seqid2taxid.map")
    for seqid, taxid in NT_TAXIDS.items():
        assert taxids[seqid] == str(taxid)
    assert UNANNOTATED_SEQUENCE not in taxids
    assert open(str(packed) + ".unannotated").read().split() == [UNANNOTATED_SEQUENCE]


def test_resume_answer_adds_nothing(tmp_path, monkeypatch):
    """The other interactive answer of a rerun keeps the library as it is: the packed file that is
    already in it must not be appended a second time."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed)
    headers = library_headers(dbname)
    assert len(headers) == 9

    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, nt_source=packed, answers=["n", "u"])
    assert library_headers(dbname) == headers


def test_multifile_in_the_genome_folder_gets_its_map_too(tmp_path, monkeypatch):
    """The other packed route: --multifile_prefix finds the archive in the genome folder, its sequences
    used to be cat'd onto the library with no map line at all."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    shutil.move(str(packed), str(genomes / "nt.fna.gz"))
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database,
                                    genomes_path=genomes, multifile_prefix="nt")

    assert len(library_headers(dbname)) == 9
    assert map_lines(out / "seqid2taxid.map") == {seqid: str(taxid) for seqid, taxid in NT_TAXIDS.items()}
    assert UNANNOTATED_SEQUENCE not in build_map(out, tmpdir)
    assert open(str(genomes / "nt.fna.gz") + ".unannotated").read().split() == [UNANNOTATED_SEQUENCE]


def test_the_header_pass_names_the_start_method_it_needs(tmp_path, monkeypatch):
    """The children of the parallel header pass inherit the database connection and the genome paths of
    the build instead of having them pickled to them, so the pass cannot rely on the default start method
    of the interpreter: python 3.14 made forkserver the default, and that pickles the target and dies on
    the sqlite connection (issue #130). The pass asks for fork itself, whatever the default is."""
    requested = []
    real_get_context = multiprocessing.get_context

    def spy(name=None):
        requested.append(name)
        return real_get_context(name)

    monkeypatch.setattr(multiprocessing, "get_context", spy)
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    out, tmpdir, dbname = run_build(tmp_path, monkeypatch, database, genomes_path=genomes)

    assert "fork" in requested
    assert len(library_headers(dbname)) == 4


@pytest.mark.parametrize("returncode", [127, 1])
def test_a_failing_program_build_is_not_reported_as_a_created_database(tmp_path, monkeypatch, returncode):
    """The status of the call of the external program decides whether there is a database: a kraken2 that
    is not installed (the shell answers 127, "not found") or a build step of its own that fails used to
    print its error on stderr while flextaxd logged that the database was created and left a database
    directory holding library.fna, prelim_map.txt and an empty inspect.txt.gz (issue #132)."""
    database = make_database(tmp_path)
    genomes, packed = copy_inputs(tmp_path)
    fakebin = tmp_path / "bin"
    fakebin.mkdir()
    failing = fakebin / "kraken2-build"
    failing.write_text("#!/bin/sh\necho 'kraken2-build: build failed' >&2\nexit {code}\n".format(code=returncode))
    failing.chmod(0o755)
    """the real programs stay reachable, only the build step of kraken2 is replaced by one that fails"""
    monkeypatch.setenv("PATH", str(fakebin) + os.pathsep + os.environ["PATH"])

    with pytest.raises(SystemExit) as stop:
        run_build(tmp_path, monkeypatch, database, genomes_path=genomes, nt_source=packed)
    assert stop.value.code == 1

    """the library step had finished: what is missing is the database the failed program was supposed to
    write, and no inspect file of a database that does not exist"""
    assert os.path.exists(str(tmp_path / "build" / "k2db" / "library" / "library.fna"))
    assert not glob.glob(str(tmp_path / "build" / "k2db" / "*.k2d"))
    assert not glob.glob(str(tmp_path / "build" / "k2db" / "inspect.txt*"))
