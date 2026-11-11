"""issue #151: DownloadGenomes imported its own package by the top level name `modules`, so the file
could not be imported under the name it really has, `flextaxd.modules.DownloadGenomes`. The command
line never noticed because custom_taxonomy_databases.py:54 puts the package folder itself on sys.path
(issue #109), but anything that imports the file by its package name - a child process of a download
run, a test, an editor - hits that failure.

Every module of the package is imported here under its own name, so a new import written the same way
fails the suite.
"""

import importlib
import pkgutil

import flextaxd


def test_every_module_of_the_package_imports_under_its_own_name():
    names = [module.name for module in pkgutil.walk_packages(flextaxd.__path__, "flextaxd.")]
    assert "flextaxd.modules.DownloadGenomes" in names
    for name in names:
        importlib.import_module(name)
