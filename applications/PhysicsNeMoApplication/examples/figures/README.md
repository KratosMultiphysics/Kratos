# Documentation figures

`make_docs_figures.py` generates the *scripted* figures of the documentation
pages under `docs/pages/Applications/PhysicsNeMo_Application/`: the
tessellation and higher-order-mode illustrations, the mesh graph and its
bistride hierarchy, the particle proximity graphs, grid sampling, the
calibration views, the halo partition, and the benchmark cost chart. Every one
of them uses the application's numpy-only bridge paths or synthetic numpy data
with matplotlib - no torch, no physicsnemo, no GPU - so the script runs in the
torch-free CI and `tests/test_docs_figures.py` executes a subset of it.

```bash
python3 make_docs_figures.py --list
python3 make_docs_figures.py                       # regenerate all into the docs tree
python3 make_docs_figures.py --only halo_partition --out /tmp/figs
```

The concept diagrams (architecture, lifecycle, data model, ...) are
hand-authored SVG files committed next to the pages; the pages reference the
SVG, which renders on the site and on GitHub without any script. A PNG twin of
every SVG is kept alongside it for viewers that cannot render SVG - both
formats are committed. After editing an SVG, regenerate its PNG with

```bash
python3 make_docs_figures.py --rasterize --only calibration_views   # PNG twins of every SVG, plus one figure
```

(`--rasterize` needs `cairosvg`; `--only` keeps the scripted figures untouched).
