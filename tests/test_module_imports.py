"""issue #151: DownloadGenomes imported its own package by the top level name `modules`, so the file
could not be imported under the name it really has, `flextaxd.modules.DownloadGenomes`. The command
line never noticed because custom_taxonomy_databases.py:54 puts the package folder itself on sys.path
(issue #109), but anything that imports the file by its package name - a child process of a download
run, a test, an editor - hits that failure.

Every module of the package is imported here under its own name, so a new import written the same way
fails the suite.
"""

import importlib
import os
import pkgutil
import subprocess
import sys

import flextaxd


def test_every_module_of_the_package_imports_under_its_own_name():
    names = [module.name for module in pkgutil.walk_packages(flextaxd.__path__, "flextaxd.")]
    assert "flextaxd.modules.DownloadGenomes" in names
    for name in names:
        importlib.import_module(name)


CHILD = """
import sys, flextaxd
sys.argv = ["flextaxd", "-db", {db!r}, "-tf", {tree!r}, "-tt", "MP", "--force"]
from flextaxd.custom_taxonomy_databases import main
main()
leaked = sorted(name for name in sys.modules if name == "modules" or name.startswith("modules."))
assert not leaked, leaked
package_folder = flextaxd.__path__[0].rstrip("/")
assert package_folder not in [p.rstrip("/") for p in sys.path], sys.path
"""


def test_the_command_does_not_make_the_package_folder_a_top_level_name(tmp_path):
    """issue #109: custom_taxonomy_databases.py put the folder of the package on sys.path, which is what
    let its own "from modules.database.DatabaseConnection import ..." lines work. That turns modules,
    database and functions into top level names for every process that imports flextaxd, and it is why a
    module importing its own package as "modules" (issue #151) went unnoticed: run the command over a real
    import and neither the name "modules" nor the folder itself may appear."""
    tree = tmp_path / "taxonomy.tsv"
    tree.write_text("RS_GCF_000005845.2\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
                    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n")
    code = CHILD.format(db=str(tmp_path / "path.fdb"), tree=str(tree))
    environment = dict(os.environ)
    repository = os.path.dirname(flextaxd.__path__[0])
    environment["PYTHONPATH"] = os.pathsep.join(
        [repository] + ([environment["PYTHONPATH"]] if environment.get("PYTHONPATH") else []))
    run = subprocess.run([sys.executable, "-c", code], cwd=str(tmp_path), env=environment,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert run.returncode == 0, run.stdout.decode()
    assert "Traceback" not in run.stdout.decode()
