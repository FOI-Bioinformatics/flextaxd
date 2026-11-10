import gzip
import os
import shutil
import sqlite3

from flextaxd.modules.ReadTaxonomyNCBI import ReadTaxonomyNCBI

TESTDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testdata", "ncbi")


def ncbi_reader(tmp_path, database):
    """ReadTaxonomyNCBI connected to a database holding the NCBI test taxonomy"""
    reader = ReadTaxonomyNCBI(taxonomy_file=os.path.join(TESTDATA, "nodes.dmp"),
                              database=str(database),
                              force_multisource=False)
    reader.parse_taxonomy()
    return reader


def genome_rows(database):
    conn = sqlite3.connect(str(database))
    return sorted(conn.execute("SELECT genome,id,reference FROM genomes").fetchall())


def annotation_nodes_missing(database):
    """The check the library step will get in T4, here as a guard for the tests: no genome
    annotation may point at an id that is not a node in the tree"""
    conn = sqlite3.connect(str(database))
    return set(row[0] for row in conn.execute(
        "SELECT DISTINCT id FROM genomes WHERE id NOT IN (SELECT id FROM nodes)"))


def write_annotation(path, header, rows):
    with gzip.open(str(path), "wt") as annotation:
        annotation.write("\t".join(header) + "\n")
        for row in rows:
            annotation.write("\t".join(row) + "\n")
    return str(path)


def write_genome(path, seqid, description="genome"):
    with open(str(path), "w") as genome_file:
        genome_file.write(">{} {}\nACGTACGT\n".format(seqid, description))
    return str(path)


def test_plain_and_gzipped_genome_files_annotate_identically(tmp_path):
    """The sequence names of the genome files are compared as bytes, so the compression sniffing
    opener has to keep the data binary: plain .fna files annotate like gzipped ones (issue #111)"""
    plain = tmp_path / "genomes_plain"
    plain.mkdir()
    for entry in os.listdir(os.path.join(TESTDATA, "genomes_plain")):
        if not entry.startswith("nt"):
            shutil.copy(os.path.join(TESTDATA, "genomes_plain", entry), str(plain / entry))
    annotation = os.path.join(TESTDATA, "gi.accession2taxid.gz")

    expected = [("GCF_000005845.2", 511145, "refseq"),
                ("GCF_000009065.1", 1841072, "refseq"),
                ("GCF_000013465.1", 367579, "refseq")]
    results = {}
    for label, folder in (("gzipped", os.path.join(TESTDATA, "genomes")), ("plain", str(plain))):
        database = tmp_path / "{}.fdb".format(label)
        reader = ncbi_reader(tmp_path, database)
        reader.parse_genomeid2taxid(folder, annotation)
        results[label] = genome_rows(database)
        assert annotation_nodes_missing(database) == set()
    assert results["gzipped"] == expected
    assert results["plain"] == expected


def test_old_annotation_file_layout_uses_the_taxid_column(tmp_path):
    """The old nucl_gb.accession2taxid layout holds the GI in the column before the taxid: the
    columns are taken from the header, the GI is never stored as a taxid (issue #110)"""
    database = tmp_path / "oldlayout.fdb"
    reader = ncbi_reader(tmp_path, database)
    reader.parse_genomeid2taxid(os.path.join(TESTDATA, "genomes"),
                               os.path.join(TESTDATA, "nucl_gb.accession2taxid.gz"))
    rows = genome_rows(database)
    assert rows == [("GCF_000005845.2", 511145, "refseq"),
                    ("GCF_000009065.1", 1841072, "refseq"),
                    ("GCF_000013465.1", 367579, "refseq")]
    stored = [str(row[1]) for row in rows]
    for gi in ("127929822", "203906411", "203906532"):
        assert gi not in stored
    assert annotation_nodes_missing(database) == set()


def test_annotation_columns_are_taken_from_the_header(tmp_path):
    """Every accession2taxid layout is identified by its header, a file without a taxid column is
    refused instead of annotated with the numbers of another column (issue #110)"""
    reader = ncbi_reader(tmp_path, tmp_path / "columns.fdb")
    assert reader.get_annotation_columns(
        [b"accession", b"accession.version", b"taxid", b"gi"]) == (1, 2)
    assert reader.get_annotation_columns(
        [b"accessions", b"accessions.versioned",
         b"accessions.nonredundant_GI_seqVersion", b"taxid"]) == (1, 3)
    assert reader.get_annotation_columns([b"gi", b"accession.version", b"taxid"]) == (1, 2)
    try:
        reader.get_annotation_columns([b"gi", b"accession.version", b"seqversion"])
        assert False, "an annotation file without a taxid column was accepted"
    except TypeError:
        pass


def test_reference_is_recorded_per_genome_file(tmp_path):
    """A genome folder that mixes assemblies with an nt archive labels every genome with the
    reference of the file it came from, the nt label does not spread to the assemblies (issue #113)"""
    genomes = tmp_path / "mixed"
    genomes.mkdir()
    write_genome(genomes / "GCF_000005845.2_ASM584v2_genomic.fna", "NC_000913.3")
    write_genome(genomes / "nt.fna", "NR_024596.1", "nt entry")
    annotation = write_annotation(tmp_path / "nucl_gb.accession2taxid.gz",
                                  ["accession", "accession.version", "taxid", "gi"],
                                  [["NC_000913", "NC_000913.3", "511145", "127929822"],
                                   ["NR_024596", "NR_024596.1", "562", "371452504"]])
    database = tmp_path / "mixed.fdb"
    reader = ncbi_reader(tmp_path, database)
    reader.parse_genomeid2taxid(str(genomes), annotation)
    assert genome_rows(database) == [("GCF_000005845.2", 511145, "refseq"),
                                     ("NR_024596.1", 562, "nt")]
    assert annotation_nodes_missing(database) == set()


def test_genome_header_without_description_still_matches(tmp_path):
    """A FASTA header may hold no description, the id is the first word of the line and not
    everything before the first space, which left the line ending in the id (issue #119)"""
    genomes = tmp_path / "headers"
    genomes.mkdir()
    with open(str(genomes / "GCF_000005845.2_ASM584v2_genomic.fna"), "w") as genome_file:
        genome_file.write(">NC_000913.3\nACGTACGT\n")
    annotation = write_annotation(tmp_path / "nucl_gb.accession2taxid.gz",
                                  ["accession", "accession.version", "taxid", "gi"],
                                  [["NC_000913", "NC_000913.3", "511145", "127929822"]])
    database = tmp_path / "headers.fdb"
    reader = ncbi_reader(tmp_path, database)
    missing = reader.parse_genomeid2taxid(str(genomes), annotation)
    assert missing == []
    assert genome_rows(database) == [("GCF_000005845.2", 511145, "refseq")]


def test_nothing_missing_writes_no_report_file(tmp_path):
    """A run that found nothing missing leaves no empty FlexTaxD.not_added behind (issue #120)"""
    assert ncbi_reader(tmp_path, tmp_path / "empty.fdb").write_missing([]) is False
    assert not os.path.exists("FlexTaxD.not_added")
