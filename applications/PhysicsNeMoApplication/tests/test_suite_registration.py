"""Guards five things that rot while the suite stays green.

The first two are static: they parse sources with `ast` and never import
the tests or the application, so they run unchanged in the torch-free CI
and cannot be disabled by a missing optional dependency. The third runs the
benchmark script, which is pure Kratos + numpy for the same reason.

**Registration.** Suites are assembled by hand, and CLAUDE.md lists the
consequence as a known gotcha - adding a `test_*.py` does not add it to the
suite. An unregistered class is invisible: the run passes, the count looks
healthy, and the tests never execute. Nothing else in the suite notices.

**Documented identifiers.** The README and the documentation pages name
modules and functions. A rename leaves the prose pointing at something that
no longer exists, and no test would fail.

**Documented import paths.** The identifier check above matches modules by
*basename*, so it survives - and therefore cannot catch - a module moving
between packages. Since python_scripts/ is a tree of packages, a wrong
`KratosMultiphysics.PhysicsNeMoApplication.<path>` is the likeliest stale
reference there is, and the `"kratos_module"`/`"python_module"` pair that
`ProjectParameters.json` uses to build a process is a path too: naming the
wrong package there fails at run time, in a user's case, not here.

**Cross-module attribute access.** Resolving an import proves the module is
importable, not that the name you then reach for is still in it - so
*splitting* a module moves names without breaking a single import statement,
and every caller keeps importing successfully and fails at the attribute.
That is not hypothetical: it is what splitting the streaming process out of
`streaming_dataset` did to its own test, which the torch-free CI would never
have run.

**The benchmark.** `benchmarks/benchmark_bridges.py` is executed by nothing
- not CMake, not this suite, not CI - while the README cites its result to
justify a roadmap decision. See the class docstring.

What the identifier check deliberately does *not* do: catch a stale claim about
project state. The roadmap once said the DoMINO predictor-corrector recipe
was still open after it had shipped; every identifier in that sentence
existed, so this check would have passed it. Prose about status has to be
read against the tree. This catches renames, which is a different failure.
"""

import ast
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest

_TESTS_DIR = Path(__file__).resolve().parent
_APP_DIR = _TESTS_DIR.parent
_SCRIPTS_DIR = _APP_DIR / "python_scripts"
_BENCHMARK = _APP_DIR / "benchmarks" / "benchmark_bridges.py"
_DOCS_DIR = (_APP_DIR.parent.parent / "docs" / "pages" / "Applications"
             / "PhysicsNeMo_Application")

_RUNNERS = ("test_PhysicsNeMoApplication.py", "test_PhysicsNeMoApplication_mpi.py")

_PACKAGE = "KratosMultiphysics.PhysicsNeMoApplication"

# tokens ending in one of these are file names, not attribute references:
# `graph_bridge.py` must not be read as module graph_bridge, attribute py
_FILE_SUFFIXES = (".py", ".md", ".json", ".mdpa", ".npz", ".pt", ".mdlus",
                  ".h5", ".so", ".txt", ".sh", ".bat", ".onnx", ".ipynb",
                  ".vtu", ".zarr", ".stl", ".csv", ".yaml", ".yml")


def _ParseFile(path):
    try:
        return ast.parse(Path(path).read_text(errors="ignore"))
    except (SyntaxError, OSError):
        return None


def _TestCaseClassesByFile():
    """{class name: file name} for every TestCase subclass outside the runners."""
    found = {}
    # rglob, not glob: the suite puts subdirectories such as mesh_bridge/ on
    # sys.path, so their test files are as much part of the suite as the top
    # level ones - and were the blind spot that let one slip through
    for path in sorted(_TESTS_DIR.rglob("test_*.py")):
        if path.name in _RUNNERS or "__pycache__" in path.parts:
            continue
        tree = _ParseFile(path)
        if tree is None:
            continue
        for node in tree.body:
            if not isinstance(node, ast.ClassDef):
                continue
            bases = [getattr(b, "attr", None) or getattr(b, "id", "")
                     for b in node.bases]
            if any("TestCase" in base for base in bases):
                found[node.name] = path.name
    return found


def _RegisteredClasses():
    """Class names actually added to a suite.

    Only names inside a `loadTestsFromTestCases([...])` call count. An import
    alone does not: importing a class and forgetting to add it is precisely
    the half-finished registration this is here to catch.
    """
    registered = set()
    for name in _RUNNERS:
        path = _TESTS_DIR / name
        if not path.is_file():
            continue
        for group in re.findall(r"loadTestsFromTestCases\(\[([\w,\s]+)\]\)",
                                path.read_text(errors="ignore")):
            registered |= {token.strip() for token in group.split(",") if token.strip()}
    return registered


def _PublicNamesByModule():
    """{module stem: {public top-level def/class names}} for the application.

    Read with `ast` rather than by importing, so a module whose imports need
    torch is still covered here.
    """
    public = {}
    for path in glob.glob(str(_SCRIPTS_DIR / "**" / "*.py"), recursive=True):
        stem = os.path.basename(path)[:-3]
        if stem == "__init__":
            continue
        tree = _ParseFile(path)
        if tree is None:
            continue
        public[stem] = {node.name for node in tree.body
                        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                             ast.ClassDef))
                        and not node.name.startswith("_")}
    return public


def _DocumentationFiles():
    files = [_APP_DIR / "README.md"]
    files += [Path(p) for p in
              glob.glob(str(_DOCS_DIR / "**" / "*.md"), recursive=True)]
    return [f for f in files if f.is_file()]


def _TopLevelNames(path):
    """Every top-level name a module binds, not only the public ones."""
    tree = _ParseFile(path)
    if tree is None:
        return set()
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names |= {target.id for target in node.targets
                      if isinstance(target, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names |= {alias.asname or alias.name.split(".")[0] for alias in node.names}
    return names


def _ResolveDottedPath(dotted):
    """Maps a path below the application package onto the real tree.

    Returns (kind, path, leftover) where kind is "package", "module" or None
    and leftover is whatever the tree could not account for - attributes of a
    module, or a path that does not exist at all.
    """
    parts = [part for part in dotted.split(".") if part]
    current, kind, path = _SCRIPTS_DIR, "package", _SCRIPTS_DIR
    for index, part in enumerate(parts):
        if (current / part).is_dir():
            current = path = current / part
            kind = "package"
            continue
        if (current / f"{part}.py").is_file():
            path = current / f"{part}.py"
            return "module", path, parts[index + 1:]
        return None, path, parts[index:]
    return kind, path, []


def _ReferenceFiles():
    """Everything that can name a module path: prose, notebooks and code."""
    patterns = ("README.md", "python_scripts/README.md", "python_scripts/**/*.py",
                "tests/**/*.py", "examples/notebooks/*.ipynb", "benchmarks/*.py")
    files = []
    for pattern in patterns:
        files += [Path(p) for p in glob.glob(str(_APP_DIR / pattern), recursive=True)]
    files += [Path(p) for p in glob.glob(str(_DOCS_DIR / "**" / "*.md"), recursive=True)]
    return [f for f in files if f.is_file() and "__pycache__" not in f.parts]


def _ApplicationSourceFiles():
    """Everything that imports the application: sources, tests, benchmarks."""
    files = []
    for pattern in ("python_scripts/**/*.py", "tests/**/*.py", "benchmarks/*.py"):
        files += [Path(p) for p in glob.glob(str(_APP_DIR / pattern), recursive=True)]
    return [f for f in files if f.is_file() and "__pycache__" not in f.parts]


def _ImportedModuleAliases(tree):
    """{local name: module file} for modules imported out of the application."""
    aliases = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.ImportFrom) and node.module
                and node.module.startswith(_PACKAGE)):
            continue
        dotted = node.module[len(_PACKAGE):].lstrip(".")
        for alias in node.names:
            kind, path, leftover = _ResolveDottedPath(
                f"{dotted}.{alias.name}" if dotted else alias.name)
            if kind == "module" and not leftover:
                aliases[alias.asname or alias.name] = path
    return aliases


class TestSuiteRegistration(KratosUnittest.TestCase):

    def test_EveryTestCaseClassIsRegisteredInASuite(self):
        defined = _TestCaseClassesByFile()
        self.assertTrue(
            defined,
            f"no TestCase classes found under {_TESTS_DIR}; this guard is not "
            "actually checking anything")

        unregistered = sorted(set(defined) - _RegisteredClasses())
        self.assertFalse(
            unregistered,
            "these TestCase classes are never added to a suite, so their tests "
            "do not run - add them to AssembleTestSuites() in "
            f"{_RUNNERS[0]} or {_RUNNERS[1]}: "
            + ", ".join(f"{defined[name]}::{name}" for name in unregistered))


class TestDocumentedIdentifiersExist(KratosUnittest.TestCase):

    def test_BacktickedModulesAndAttributesResolve(self):
        public = _PublicNamesByModule()
        self.assertTrue(public, f"no modules found under {_SCRIPTS_DIR}")
        documents = _DocumentationFiles()
        self.assertTrue(documents, "no documentation files found to check")

        broken = []
        module_hits = attribute_hits = 0
        for document in documents:
            tokens = set(re.findall(r"`([^`\n]+)`",
                                    document.read_text(errors="ignore")))
            for token in tokens:
                name = token.strip().replace("()", "")
                if name.endswith(_FILE_SUFFIXES):
                    continue
                if name in public:
                    module_hits += 1
                    continue
                dotted = re.fullmatch(r"([a-z_][a-z0-9_]*)\.([A-Za-z_]\w*)", name)
                # only our own modules; external APIs are not ours to verify
                if dotted and dotted.group(1) in public:
                    attribute_hits += 1
                    if dotted.group(2) not in public[dotted.group(1)]:
                        broken.append(f"{document.name}: `{token}`")

        self.assertFalse(
            broken,
            "documentation names application code that does not exist "
            "(renamed or removed?): " + ", ".join(sorted(set(broken))))
        self.assertTrue(
            module_hits or attribute_hits,
            "no documented references were resolved at all, so this guard "
            "passed without checking anything")


class TestDocumentedImportPathsResolve(KratosUnittest.TestCase):
    """Every written-down module path must exist in the tree.

    `TestDocumentedIdentifiersExist` matches modules by basename, so moving
    one between packages leaves it green while every documented import breaks.
    This resolves the *path*, which is the thing the reorganisation can get
    wrong.
    """

    def test_EveryWrittenModulePathExists(self):
        pattern = re.compile(rf"{re.escape(_PACKAGE)}((?:\.\w+)*)")

        broken, resolved = [], 0
        for document in _ReferenceFiles():
            for match in pattern.finditer(document.read_text(errors="ignore")):
                dotted = match.group(1).lstrip(".")
                kind, _, leftover = _ResolveDottedPath(dotted)
                resolved += 1
                if kind is None:
                    broken.append(f"{document.name}: {match.group(0)} "
                                  f"(no '{leftover[0]}' there)")
                elif kind == "package" and leftover:
                    broken.append(f"{document.name}: {match.group(0)} "
                                  f"(no module '{leftover[0]}' in that package)")

        self.assertFalse(
            broken,
            "these module paths do not exist - a package moved without the "
            "reference following it: " + ", ".join(sorted(set(broken))))
        self.assertTrue(resolved > 100,
                        f"only {resolved} paths seen; this guard is not "
                        "actually reading the tree")

    def test_EveryWrittenAttributePathExists(self):
        """The tail of a path, when there is one, must be a real name."""
        pattern = re.compile(rf"{re.escape(_PACKAGE)}((?:\.\w+)+)")

        broken = []
        for document in _ReferenceFiles():
            for match in pattern.finditer(document.read_text(errors="ignore")):
                kind, path, leftover = _ResolveDottedPath(match.group(1).lstrip("."))
                if kind != "module" or not leftover:
                    continue
                if leftover[0] not in _TopLevelNames(path):
                    broken.append(f"{document.name}: {match.group(0)} "
                                  f"('{path.name}' defines no '{leftover[0]}')")

        self.assertFalse(
            broken,
            "these paths name something their module does not define: "
            + ", ".join(sorted(set(broken))))

    def test_EveryDocumentedProcessFactoryIsBuildable(self):
        """`kratos_module` + `python_module` is a path too, and a fragile one.

        Kratos builds a process by importing `<kratos_module>.<python_module>`
        and calling its `Factory`. Getting the package wrong there fails in a
        user's case, at run time, with nothing here to catch it - so the pair
        is resolved and the target checked for the `Factory` that makes it
        usable at all.
        """
        pair = re.compile(
            r'"python_module"\s*:\s*"(\w+)"(?:(?!"python_module").)*?'
            r'"kratos_module"\s*:\s*"([\w.]+)"'
            r'|"kratos_module"\s*:\s*"([\w.]+)"(?:(?!"kratos_module").)*?'
            r'"python_module"\s*:\s*"(\w+)"', re.DOTALL)

        broken, checked = [], 0
        for document in _ReferenceFiles():
            for match in pair.finditer(document.read_text(errors="ignore")):
                module = match.group(1) or match.group(4)
                package = match.group(2) or match.group(3)
                if not package.startswith(_PACKAGE):
                    continue          # core and other applications are not ours
                checked += 1
                dotted = f"{package[len(_PACKAGE):].lstrip('.')}.{module}".lstrip(".")
                kind, path, _ = _ResolveDottedPath(dotted)
                if kind != "module":
                    broken.append(f"{document.name}: \"{package}\" has no "
                                  f"\"{module}\"")
                elif "Factory" not in _TopLevelNames(path):
                    broken.append(f"{document.name}: {path.name} defines no "
                                  "Factory, so it cannot be built from "
                                  "ProjectParameters")

        self.assertFalse(
            broken,
            "these documented process declarations cannot be instantiated: "
            + ", ".join(sorted(set(broken))))
        self.assertTrue(checked, "no process declarations found to check")


class TestCrossModuleAttributesExist(KratosUnittest.TestCase):
    """`from ... import mod` then `mod.Name` - does mod still define Name?

    Import resolution stops one step short of this, and that step is exactly
    where splitting a module hurts: the names move, every import keeps
    working, and the failure only appears when something calls the attribute.
    In a suite whose ML tests skip without torch, that can be never.
    """

    def test_EveryAttributeOfAnImportedModuleIsDefined(self):
        broken, checked = [], 0
        for path in _ApplicationSourceFiles():
            tree = _ParseFile(path)
            if tree is None:
                continue
            aliases = _ImportedModuleAliases(tree)
            if not aliases:
                continue
            for node in ast.walk(tree):
                if not (isinstance(node, ast.Attribute)
                        and isinstance(node.value, ast.Name)):
                    continue
                target = aliases.get(node.value.id)
                if target is None:
                    continue
                checked += 1
                if node.attr not in _TopLevelNames(target):
                    broken.append(
                        f"{path.name}:{node.lineno}: {node.value.id}.{node.attr} "
                        f"- {target.name} defines no \"{node.attr}\"")

        self.assertFalse(
            broken,
            "these call a name their module does not define - moved to "
            "another module? " + ", ".join(sorted(set(broken))))
        self.assertTrue(checked > 100,
                        f"only {checked} attribute accesses seen; this guard "
                        "is not actually reading the sources")


class TestBenchmarkStillRuns(KratosUnittest.TestCase):
    """The benchmark is a script nothing else executes.

    Not CMake, not this suite, not CI - and the README cites its result to
    justify a roadmap decision. That is how its headline number drifted from
    what it actually measures without anyone noticing, so at minimum the
    script itself must keep running against the bridges it times.

    Deliberately asserts nothing about *timings*. A wall-clock assertion in
    a unit suite is flaky by construction; this checks the API, and the
    numbers stay a thing a human runs and reads.
    """

    def test_TheBridgeBenchmarkExecutes(self):
        self.assertTrue(_BENCHMARK.is_file(), f"{_BENCHMARK} is missing")

        completed = subprocess.run(
            [sys.executable, str(_BENCHMARK),
             "--divisions", "4", "--grid", "8", "--repeat", "1"],
            capture_output=True, text=True, timeout=600)

        self.assertEqual(
            completed.returncode, 0,
            f"{_BENCHMARK.name} failed:\n{completed.stdout[-2000:]}\n"
            f"{completed.stderr[-2000:]}")
        self.assertIn("us/entity", completed.stdout,
                      "the benchmark ran but printed no results table")


if __name__ == "__main__":
    KratosUnittest.main()
