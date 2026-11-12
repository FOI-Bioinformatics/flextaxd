"""T4: the library of a build, counted before the program reads it.

kraken2, krakenuniq and ganon take the taxid of a sequence from the library header or from a line of
the taxid map keyed by the id the library holds, and a sequence that reaches neither is left out of
the program database without the program saying so. The numbers below are the ones a build reports
before it calls the program: how many sequences and map lines it assembled, how many of its sequences
the program can give no taxid to, and how many map lines name an annotated genome the run had no
sequence file of its own. That last number is what the disk report of issue #131 makes trustworthy:
the sequence level annotations of a packed file have no file of their own by design and must not show
up as genomes without a file.

These tests are made of hand-written files: only the kraken route has a place where the counting runs
today, since which programs a build supports is an open question (#150).
"""

import os

from flextaxd.modules.database.DatabaseConnection import ModifyFunctions
from flextaxd.modules.library_check import (library_numbers, library_report, map_key, sequence_key,
                                            store_library_numbers)


def write(path, lines):
    with open(str(path), "w") as handle:
        for line in lines:
            handle.write(line + "\n")


def kraken_files(tmp_path):
    """A library and a prelim map in the shape kraken2 is given them: four sequences (two with the
    kraken taxid tag, one annotated through the map, one with no annotation at all) and five map lines,
    two of which name annotated genomes that have no sequence file"""
    library = tmp_path / "library.fna"
    write(library, [">BA000036.1|kraken:taxid|1280  ", "ACGTACGT",
                    ">NC_000913.3 Escherichia coli str. K-12, complete genome", "ACGTACGT",
                    ">CP000250.1 Escherichia coli O157:H7, complete genome", "ACGTACGT",
                    ">LC_test_1 sequence with no annotation anywhere", "ACGTACGT"])
    map_file = tmp_path / "prelim_map.txt"
    write(map_file, ["TAXID\tBA000036.1|kraken:taxid|1280\t1280",
                     "TAXID\tNC_000913.3\t511145",
                     "TAXID\tCP000250.1\t562",
                     "TAXID\tGCF_000001415.1\t1423",
                     "TAXID\tGCF_000009605.1\t28901"])
    """the database annotated six genomes, two of them have no file on disk; the four sequences the
    packed file brought in are annotated on sequence level"""
    annotated = ["BA000036.1", "GCF_000005845.2", "GCF_000013465.1",
                 "GCF_000001415.1", "GCF_000009605.1", "NC_000913.3", "CP000250.1"]
    files = ["BA000036.1", "GCF_000005845.2", "GCF_000013465.1"]
    sequence_level = ["NC_000913.3", "CP000250.1"]
    return library, map_file, annotated, files, sequence_level


def test_the_library_and_its_map_are_counted_as_the_program_reads_them(tmp_path):
    """Four sequences, five map lines, one sequence the program can give no taxid to, and the two map
    lines naming a genome the run had no file of."""
    library, map_file, annotated, files, sequence_level = kraken_files(tmp_path)

    numbers = library_numbers([str(library)], str(map_file),
                              annotated_genomes=annotated,
                              genome_files=files,
                              sequence_level_annotations=sequence_level,
                              map_key_column=1)

    assert numbers["sequences"] == 4
    assert numbers["map_lines"] == 5
    assert numbers["without_taxid"] == 1
    assert numbers["map_without_genome"] == 2


def test_a_genome_sequence_needs_only_a_tag_or_a_map_line(tmp_path):
    """The sequence with the |kraken:taxid| tag and the plain accession that has a map line are both
    seen by kraken2; only the sequence with neither counts as without a taxid. Dropping the map line of
    the plain accession is what makes it one."""
    library, map_file, annotated, files, sequence_level = kraken_files(tmp_path)
    write(map_file, [line for line in open(str(map_file))
                     if not line.startswith("TAXID\tCP000250.1")])

    numbers = library_numbers([str(library)], str(map_file),
                              annotated_genomes=annotated,
                              genome_files=files,
                              sequence_level_annotations=sequence_level,
                              map_key_column=1)

    assert numbers["sequences"] == 4
    assert numbers["without_taxid"] == 2, "the accession that lost its map line is now without a taxid"


def test_sequence_level_annotations_are_not_genomes_without_a_file(tmp_path):
    """Issue #131 is what makes the fourth number mean what it says. The sequences of a packed file are
    annotated on sequence level and their taxids are map lines; their accession is a genome name in the
    genomes table, so counting the annotations without a file of their own would name four map lines
    that are perfectly right."""
    library, map_file, annotated, files, sequence_level = kraken_files(tmp_path)

    without = library_numbers([str(library)], str(map_file),
                              annotated_genomes=annotated,
                              genome_files=files,
                              sequence_level_annotations=(),
                              map_key_column=1)
    assert without["map_without_genome"] == 4, "without the sequence level names the count is wrong by four"

    with_sequence_level = library_numbers([str(library)], str(map_file),
                                          annotated_genomes=annotated,
                                          genome_files=files,
                                          sequence_level_annotations=sequence_level,
                                          map_key_column=1)
    assert with_sequence_level["map_without_genome"] == 2


def test_a_map_keyed_by_its_first_column_is_counted_too(tmp_path):
    """The counting does not assume the kraken2 field order: a seq-info map of three columns keyed by
    the first one (the shape ganon and centrifuge are given) is read the same way, and a sequence the
    map does not name is the one the program will leave out. Written from hand-made files: no build of
    such a program is tested here, which programs the build supports is an open question (#150)."""
    per_process = tmp_path / ".tmp0.fasta"
    packed_library = tmp_path / ".tmppacked.fasta"
    write(per_process, [">NC_000913.3_GCF_000005845.2", "ACGTACGT",
                        ">BA000036.1_GCA_000016685.1", "ACGTACGT"])
    write(packed_library, [">CP000250.1", "ACGTACGT",
                           ">LC_test_1", "ACGTACGT"])
    map_file = tmp_path / "seqid2taxid.map"
    write(map_file, ["NC_000913.3_GCF_000005845.2\t8\t562",
                     "BA000036.1_GCA_000016685.1\t8\t1280",
                     "CP000250.1\t8\t562"])

    numbers = library_numbers([str(per_process), str(packed_library)], str(map_file),
                              annotated_genomes=["GCF_000005845.2", "GCA_000016685.1", "CP000250.1"],
                              genome_files=["GCF_000005845.2", "GCA_000016685.1"],
                              sequence_level_annotations=["CP000250.1"],
                              map_key_column=0)

    assert numbers["sequences"] == 4
    assert numbers["map_lines"] == 3
    assert numbers["without_taxid"] == 1, "ganon has no taxid for LC_test_1 and leaves it out"
    assert numbers["map_without_genome"] == 0


def test_the_numbers_are_stored_once_and_printed_once(tmp_path):
    """The report of the flextaxd database keeps one line about the library, next to the version and the
    date of the build, and a second build replaces it instead of adding to it. A database that never
    built anything has no such line."""
    database = str(tmp_path / "build.fdb")
    stored = ModifyFunctions(database=database)
    assert library_report(stored) is None, "a database that never built a library has no library line"

    store_library_numbers(stored, "kraken2", str(tmp_path / "k2db"),
                          {"sequences": 9, "map_lines": 18, "without_taxid": 1, "map_without_genome": 2})
    assert stored.get_meta("library_sequences") == "9"
    assert stored.get_meta("library_program") == "kraken2"
    line = library_report(stored)
    assert "9 sequences, 18 map lines, 1 without a taxid, 2 map lines for a genome without sequence file" in line
    report = "\n".join(stored.report_lines())
    assert report.count("Library of the last build") == 1

    store_library_numbers(stored, "ganon", str(tmp_path / "ganondb"),
                          {"sequences": 8, "map_lines": 8, "without_taxid": 0, "map_without_genome": 0})
    report = "\n".join(stored.report_lines())
    assert report.count("Library of the last build") == 1, "the numbers of the last build replace the earlier ones"
    assert "Library of the last build (ganon," in report
    assert "8 sequences, 8 map lines, 0 without a taxid" in report


def test_the_id_a_map_line_and_a_library_header_agree_on(tmp_path):
    """The map is keyed by the id the library holds: the kraken tag belongs to the header, not to the
    key, and a ganon key is the sequence id with its genome behind it."""
    assert sequence_key(">NC_000913.3|kraken:taxid|562  ") == ("NC_000913.3", True)
    assert sequence_key(">NC_000913.3 Escherichia coli, complete genome") == ("NC_000913.3", False)
    assert map_key("TAXID\tBA000036.1|kraken:taxid|1280\t1280", 1) == "BA000036.1"
    assert map_key("NC_000913.3_GCF_000005845.2\t8\t562", 0) == "NC_000913.3_GCF_000005845.2"
    assert map_key("empty\t", 1) is None
