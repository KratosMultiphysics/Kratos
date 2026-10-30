"""Every relative link and image reference in the documentation must resolve.

The documentation tree now carries some forty figures and a few hundred
cross-links between pages. A moved image or a renamed page leaves a broken
reference that no other test sees - the identifier guard checks *names*, not
*paths on disk*. This is torch-free and pure text.
"""

import glob
import re
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest

_TESTS_DIR = Path(__file__).resolve().parent
_APP_DIR = _TESTS_DIR.parent
_DOCS_DIR = (_APP_DIR.parent.parent / "docs" / "pages" / "Applications"
             / "PhysicsNeMo_Application")

_MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
_HTML_SRC = re.compile(r"""<(?:img|a)\s[^>]*?(?:src|href)=["']([^"']+)["']""")


def _References(text):
    return set(_MARKDOWN_LINK.findall(text)) | set(_HTML_SRC.findall(text))


def _IsExternal(reference):
    return reference.startswith(("http://", "https://", "mailto:", "#"))


def _Resolve(page, reference):
    target = reference.split("#", 1)[0]
    if not target:
        return page  # a pure anchor
    path = (page.parent / target).resolve()
    if path.suffix == ".html":
        path = path.with_suffix(".md")
    return path


class TestDocumentationLinks(KratosUnittest.TestCase):

    def _Pages(self):
        pages = [Path(p) for p in glob.glob(str(_DOCS_DIR / "**" / "*.md"), recursive=True)]
        pages.append(_APP_DIR / "README.md")
        pages.append(_APP_DIR / "examples" / "notebooks" / "README.md")
        return [p for p in pages if p.is_file()]

    def test_EveryRelativeReferenceResolves(self):
        pages = self._Pages()
        self.assertGreater(len(pages), 40, "the documentation tree looks empty")
        broken, checked = [], 0
        for page in pages:
            for reference in _References(page.read_text(errors="ignore")):
                if _IsExternal(reference):
                    continue
                checked += 1
                target = _Resolve(page, reference)
                if not target.exists():
                    broken.append(f"{page.relative_to(_APP_DIR.parent.parent)}: {reference}")
        self.assertFalse(broken, "broken relative references:\n  " + "\n  ".join(sorted(broken)))
        self.assertGreater(checked, 100, f"only {checked} references seen; this guard is not reading the pages")

    def test_EveryImageFileIsReferenced(self):
        """An image nobody references is dead weight in the repository.

        Every SVG diagram has a PNG twin next to it (both formats are kept on
        purpose); the twin counts as referenced through its SVG.
        """
        referenced = set()
        for page in self._Pages():
            for reference in _References(page.read_text(errors="ignore")):
                if not _IsExternal(reference):
                    referenced.add(_Resolve(page, reference))
        images = [Path(p) for p in glob.glob(str(_DOCS_DIR / "**" / "images" / "*.*"), recursive=True)]
        self.assertTrue(images, "no images found under the documentation tree")
        def is_referenced(path):
            if path.resolve() in referenced:
                return True
            # the PNG render kept next to an SVG diagram the pages reference
            return path.suffix == ".png" and path.with_suffix(".svg").resolve() in referenced
        orphans = [str(p.relative_to(_DOCS_DIR)) for p in images if not is_referenced(p)]
        self.assertFalse(orphans, "images no page references: " + ", ".join(sorted(orphans)))

    def test_FrontMatterHasNoColonsInValues(self):
        """docs/process_pages.py splits front-matter lines on ':' and keeps the first piece."""
        offending = []
        for page in self._Pages():
            if _DOCS_DIR not in page.parents:
                continue
            lines = page.read_text(errors="ignore").splitlines()
            if not lines or lines[0].strip() != "---":
                offending.append(f"{page.name}: no front matter")
                continue
            for line in lines[1:]:
                if line.strip() == "---":
                    break
                key, _, value = line.partition(":")
                if ":" in value:
                    offending.append(f"{page.name}: '{key.strip()}' value contains a colon")
                if key.strip() == "summary" and not value.strip():
                    offending.append(f"{page.name}: empty summary")
        self.assertFalse(offending, "\n  ".join([""] + sorted(offending)))


if __name__ == "__main__":
    KratosUnittest.main()
