"""the name of the lineage format, QIIME -> MP, as the publication has it.

The format is the semicolon separated lineage one genome per row that a GTDB taxonomy.tsv writes
(`d__Bacteria;p__...;g__...;s__...`). It was read by ReadTaxonomyQIIME behind `-tt QIIME`, which names a
program rather than the layout, so the reader is ReadTaxonomyMP and the option is `-tt MP`. The file and the
reader are the same code, so the two things that have to hold are that MP reads such a file and that a script
still passing QIIME is told what to write instead of being shown a list of choices.
"""

import os
import sys

import pytest

from flextaxd import custom_taxonomy_databases

LINEAGE_FILE = (
    "RS_GCF_000005845.2\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;"
    "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
    "GB_GCA_000016685.1\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;"
    "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
)


def build(tmp_path, monkeypatch, taxonomy_type):
    """the import step of the command line, run with the taxonomy type a user would write"""
    database = str(tmp_path / "mp.fdb")
    tree_file = tmp_path / "taxonomy.tsv"
    tree_file.write_text(LINEAGE_FILE)
    monkeypatch.setattr(sys, "argv", ["flextaxd", "-db", database, "-tf", str(tree_file),
                                      "-tt", taxonomy_type, "--force"])
    try:
        custom_taxonomy_databases.main()
    except SystemExit as stop:
        assert stop.code in (None, 0), stop.code
    return database


def test_an_mp_lineage_file_imports_under_the_new_name(tmp_path, monkeypatch):
    from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

    database = build(tmp_path, monkeypatch, "MP")
    conn = ModifyFunctions(database=database)
    names = {name for _, name in conn.select("SELECT id,name FROM nodes")}
    assert "Escherichia coli" in names and "Enterobacteriaceae" in names
    assert conn.select("SELECT COUNT(*) FROM genomes")[0][0] == 2


def test_GTDB_still_maps_to_the_lineage_reader(tmp_path, monkeypatch):
    """GTDB was never a reader of its own, it has always been the GTDB file of this format"""
    from flextaxd.modules.database.DatabaseConnection import ModifyFunctions

    database = build(tmp_path, monkeypatch, "GTDB")
    names = {name for _, name in ModifyFunctions(database=database).select("SELECT id,name FROM nodes")}
    assert "Escherichia coli" in names


def test_the_old_name_is_refused_naming_the_new_one(tmp_path, monkeypatch, capsys):
    database = str(tmp_path / "old.fdb")
    tree_file = tmp_path / "taxonomy.tsv"
    tree_file.write_text(LINEAGE_FILE)
    monkeypatch.setattr(sys, "argv", ["flextaxd", "-db", database, "-tf", str(tree_file), "-tt", "QIIME"])
    with pytest.raises(SystemExit) as stop:
        custom_taxonomy_databases.main()
    assert stop.value.code == 2
    assert "QIIME is gone" in capsys.readouterr().err
    assert not os.path.exists(database), "the run went on and built a database anyway"


def test_the_help_lists_the_new_name_and_not_the_old_one(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["flextaxd", "--help"])
    with pytest.raises(SystemExit):
        custom_taxonomy_databases.main()
    help_text = capsys.readouterr().out
    assert "MP" in help_text
    assert "QIIME" not in help_text
