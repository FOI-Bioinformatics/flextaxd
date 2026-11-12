"""issue #74: the two-level model has to be in the manual, not only in the wiki

The question that started this was "how do I merge two kraken databases", which is a question about the
wrong level: flextaxd merges taxonomy databases (`-db`, the sqlite file it reads back), while the
classifier databases built from them (`--db_name`, kraken2 and friends) are output only. Both commands now
carry the same four-command tour in their --help, and the README states the model and the merge example,
so this test pins the three things a user has to be able to find: the tour, the sentence that a built
classifier database cannot be read back, and the words taxonomy database / classifier database used for the
two levels.
"""

import os
import sys

import pytest

from flextaxd import create_databases, custom_taxonomy_databases

README = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "README.md")


def help_text(module, monkeypatch, capsys, command):
    monkeypatch.setattr(sys, "argv", [command, "--help"])
    with pytest.raises(SystemExit):
        module.main()
    return capsys.readouterr().out


@pytest.mark.parametrize("module,command", [(custom_taxonomy_databases, "flextaxd"),
                                            (create_databases, "flextaxd-create")])
def test_the_help_of_both_commands_tours_the_four_steps(module, command, monkeypatch, capsys):
    text = help_text(module, monkeypatch, capsys, command)
    for step in ("1 import", "2 modify", "3 report", "4 build"):
        assert step in text, "{c} help has no {s} step".format(c=command, s=step)
    assert "cannot be read back" in text and "what can be imported are its inputs" in text, text
    assert "taxonomy database" in text and "--db_name" in text, text
    assert "flextaxd-create" in text, command


def test_the_two_commands_answer_with_the_same_model(monkeypatch, capsys):
    """one text, not two that drift apart: the wording of the tour is identical in both commands"""
    first = help_text(custom_taxonomy_databases, monkeypatch, capsys, "flextaxd")
    second = help_text(create_databases, monkeypatch, capsys, "flextaxd-create")
    assert (first.split("\n\n")[-1].strip()) == (second.split("\n\n")[-1].strip())


def test_the_readme_states_the_two_levels_and_the_merge():
    with open(README, encoding="utf-8") as manual:
        text = manual.read()
    assert "## The two levels: taxonomy database and classifier database" in text
    assert "cannot be imported into flextaxd" in text
    assert "flextaxd -db life.ftd -md ncbi.ftd --parent root -gt genome2node.tsv" in text
    assert "flextaxd-create -db life.ftd --db_name kraken_run" in text
