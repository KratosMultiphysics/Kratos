"""Executes the shipped example notebooks.

`examples/notebooks` is the application's worked documentation: twenty
notebooks calling the same public helpers the processes call. Nothing
executed them until this file existed, so a renamed argument or a changed
return type broke them silently - a static import check is not enough,
because every one of those names still resolves; it is the call signatures
and the behaviour that rot.

Each notebook runs in a throwaway copy of the tree rather than in place,
for two reasons. The notebooks write artifacts next to themselves
(`output/`, `*.pt`, `*.png`), which would dirty the working tree on every
run; and notebooks 07, 16, 17 and 19 load their cases from
`tests/kratos_solver_cases` - 19 also from `tests/adjoint_cases` - so the
copy has to preserve the relative layout, not just the notebook file.

The notebooks are slow by test standards (~5 minutes for all twenty) and
are registered on the validation suite only.
"""

import os
import shutil
import tempfile
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest

_TESTS_DIR = Path(__file__).resolve().parent
_NOTEBOOK_DIR = _TESTS_DIR.parent / "examples" / "notebooks"
_CASE_DIR = _TESTS_DIR / "kratos_solver_cases"
_ADJOINT_CASE_DIR = _TESTS_DIR / "adjoint_cases"

_CELL_TIMEOUT = 900


def _TryImportNbclient():
    try:
        import nbformat
        from nbclient import NotebookClient
        return nbformat, NotebookClient
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication's notebook tests require nbclient and nbformat, "
            "which could not be imported. Install them with e.g. "
            "'pip install nbclient nbformat'.") from e


def _HaveNbclient():
    try:
        _TryImportNbclient()
        return True
    except ImportError:
        return False


def _HaveNotebookDependencies():
    """The notebooks themselves need torch, physicsnemo and matplotlib.

    Checked in-process even though the notebook runs in a kernel
    subprocess: the subprocess inherits this interpreter's environment, so
    an import that fails here fails there too - only later, and as an
    error rather than a skip.
    """
    try:
        import matplotlib  # noqa: F401
        import physicsnemo  # noqa: F401
        import torch        # noqa: F401
        return True
    except ImportError:
        return False


def _ListNotebooks():
    if not _NOTEBOOK_DIR.is_dir():
        return []
    return sorted(_NOTEBOOK_DIR.glob("*.ipynb"))


def _RunNotebook(notebook_path):
    """Executes one notebook in a temporary copy of the tree.

    Raises whatever the notebook raised; a clean return means every code
    cell ran without an exception.
    """
    nbformat, NotebookClient = _TryImportNbclient()

    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        run_dir = root / "examples" / "notebooks"
        run_dir.mkdir(parents=True)
        # 16 and 17 reach out to ../../tests/kratos_solver_cases, and 19
        # additionally to ../../tests/adjoint_cases
        shutil.copytree(_CASE_DIR, root / "tests" / "kratos_solver_cases")
        shutil.copytree(_ADJOINT_CASE_DIR, root / "tests" / "adjoint_cases")
        shutil.copy(notebook_path, run_dir / notebook_path.name)

        # the plotting notebooks must not need a display
        previous_backend = os.environ.get("MPLBACKEND")
        os.environ["MPLBACKEND"] = "Agg"
        try:
            notebook = nbformat.read(
                str(run_dir / notebook_path.name), as_version=4)
            NotebookClient(
                notebook,
                timeout=_CELL_TIMEOUT,
                kernel_name="python3",
                resources={"metadata": {"path": str(run_dir)}}).execute()
        finally:
            if previous_backend is None:
                os.environ.pop("MPLBACKEND", None)
            else:
                os.environ["MPLBACKEND"] = previous_backend


@KratosUnittest.skipUnless(
    _HaveNbclient(), "requires nbclient ('pip install nbclient nbformat')")
@KratosUnittest.skipUnless(
    _HaveNotebookDependencies(),
    "the example notebooks require torch, physicsnemo and matplotlib")
class TestNotebooks(KratosUnittest.TestCase):
    """One generated test method per notebook, so a failure names it."""

    def test_AllNotebooksAreCovered(self):
        """Guards the generation itself.

        The methods below are attached at import time from a glob. A wrong
        path would attach nothing at all and this whole file would pass
        while testing no notebook, which is the one failure mode the
        generated tests cannot report.
        """
        notebooks = _ListNotebooks()
        self.assertTrue(
            notebooks,
            f"no notebooks found in {_NOTEBOOK_DIR}; the example notebooks are "
            "not being tested at all")

        generated = {name[len("test_Notebook_"):] for name in dir(self)
                     if name.startswith("test_Notebook_")}
        self.assertEqual(generated, {path.stem for path in notebooks})


def _MakeNotebookTest(notebook_path):
    def test(self):
        _RunNotebook(notebook_path)
    test.__name__ = f"test_Notebook_{notebook_path.stem}"
    test.__doc__ = f"Executes examples/notebooks/{notebook_path.name}."
    return test


for _path in _ListNotebooks():
    _test = _MakeNotebookTest(_path)
    setattr(TestNotebooks, _test.__name__, _test)


if __name__ == "__main__":
    KratosUnittest.main()
