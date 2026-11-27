#!/usr/bin/env python3
"""Headless runner for the generator notebooks.

Executes the code cells of one or more ``.ipynb`` files in order, in a single namespace, with
the working directory set to the folder of the notebook (as a Jupyter kernel would), so that
the condition files can be regenerated without a Jupyter installation:

    python3 run_notebook.py ALM_frictional_mortar_condition/ALM_frictional_mortar_condition.ipynb

When ``nbclient`` and ``nbformat`` are available (``pip install nbclient nbformat ipykernel``) the
notebook is executed through them instead and the outputs are stored back into the file
(``--no-store`` disables the storing). Only the Python standard library is needed otherwise.
IPython magics (``%``, ``!``) are skipped in the fallback mode.
"""
import json
import os
import re
import sys
import traceback


def _ExecuteWithNbclient(path, store):
    import nbformat
    from nbclient import NotebookClient
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(notebook, timeout=None, kernel_name=notebook.metadata.get("kernelspec", {}).get("name", "python3"), resources={"metadata": {"path": os.path.dirname(path)}})
    client.execute()
    if store:
        nbformat.write(notebook, path)


def _ExecuteWithExec(path):
    with open(path, encoding="utf-8") as notebook_file:
        notebook = json.load(notebook_file)
    if notebook.get("nbformat") != 4:
        raise RuntimeError("{}: unsupported nbformat {}".format(path, notebook.get("nbformat")))
    notebook_dir = os.path.dirname(path)
    os.chdir(notebook_dir)
    sys.path.insert(0, notebook_dir)
    namespace = {"__name__": "__main__", "__builtins__": __builtins__}
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = cell["source"]
        source = "".join(source) if isinstance(source, list) else source
        source = "\n".join(line for line in source.splitlines() if not re.match(r"^\s*[%!]", line))
        try:
            exec(compile(source, "{} [cell {}]".format(path, index), "exec"), namespace)
        except SystemExit:
            raise
        except Exception:
            traceback.print_exc()
            raise SystemExit("Notebook {}: code cell {} failed".format(path, index))


def Run(path, store=True, force_exec=False):
    path = os.path.abspath(path)
    if not force_exec:
        try:
            import nbclient  # noqa: F401
            import nbformat  # noqa: F401
            _ExecuteWithNbclient(path, store)
            return
        except ImportError:
            pass
    _ExecuteWithExec(path)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("notebooks", nargs="+")
    parser.add_argument("--no-store", action="store_true", help="Do not write the outputs back into the notebook (nbclient mode)")
    parser.add_argument("--plain", action="store_true", help="Force the standard-library exec mode even if nbclient is installed")
    args = parser.parse_args(argv)
    cwd = os.getcwd()
    for notebook in args.notebooks:
        os.chdir(cwd)
        Run(notebook, store=not args.no_store, force_exec=args.plain)


if __name__ == "__main__":
    main()
