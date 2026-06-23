# Kratos MDPA Preview (VS Code extension)

Preview, organize, and manage Kratos Multiphysics `.mdpa` model-part files
directly in VS Code: a 3D mesh viewer with a navigable **ModelPart /
SubModelPart outline** whose entries are **toggleable layers**.

It is fully self-contained — a pure-TypeScript MDPA parser feeds a
[VTK.js](https://kitware.github.io/vtk-js/) viewer running in a webview. **No
Python or compiled Kratos is required.**

## Features

- **3D preview** of nodes, elements, conditions, and geometries. Volume
  elements (tet/hex/wedge/pyramid) are shown as their boundary surface;
  quadratic elements are approximated by their corner nodes.
- **Outline tree** of the entity blocks and the full SubModelPart hierarchy,
  with per-row visibility checkboxes (activate/deactivate a layer) and
  click-to-frame.
- **SubModelParts as layers** — each SubModelPart is an independently toggleable
  overlay so you can isolate inlets/outlets/boundaries.
- **Stats panel**: node/element/condition/geometry counts, bounding box,
  detected 2D/3D, and any element type names that could not be mapped.
- **Editor integration**: `mdpa` language id with `//` comments, `Begin`/`End`
  folding, and syntax highlighting. The raw text editor stays the default; open
  the preview from the editor-title button, the explorer context menu, or the
  **Open MDPA Preview** command.

## Develop

```bash
npm install
npm run compile      # bundle extension (dist/) and webview (media/) via esbuild
npm run watch        # rebuild on change
npm test             # parser unit tests (node:test) against repo fixtures
npm run typecheck    # tsc --noEmit
```

Press **F5** in VS Code to launch an Extension Development Host, then open any
`.mdpa` file (e.g. those under `applications/*/tests/`).

## Layout

| Path | Purpose |
|------|---------|
| `src/extension.ts` | Activation, command + custom-editor registration |
| `src/mdpaEditorProvider.ts` | Custom editor: parses the document, hosts the webview |
| `src/parser/` | `mdpaParser`, `geometryMap` (Kratos name → VTK cell), `types` |
| `webview/` | `main.ts` (VTK scene), `meshBuilder.ts`, `outline.ts`, `style.css` |
| `syntaxes/` | TextMate grammar for highlighting |

The Kratos name → VTK cell-type table mirrors the core
`kratos/input_output/vtk_definition.cpp` and `kratos/sources/kratos_application.cpp`.

> This is build tooling only and is intentionally outside the Kratos C++/CMake
> build.
