---
name: geo-format
description: "Use when formatting C++ or JSON files of GeoMechanicsApplication, before you report a change as done. Covers kp format, clang-format with the GeoMechanics settings, and the JSON formatting rules."
---

# Format C++ and JSON files of GeoMechanicsApplication

Paths below are relative to `applications/GeoMechanicsApplication/`. The formatting settings
are in `.clang-format`. Read them there; this skill does not repeat them.

## With a `kp` environment

When a `kp` environment is available, run `kp format`. It formats the C++ files and the new
JSON files.

## Without a `kp` environment

### C++

- Format changed C++ code with `clang-format` using
  `applications/GeoMechanicsApplication/.clang-format`. Prefer formatting only the changed lines
  or ranges to avoid unrelated whole-file changes.
- To format only the changed lines, run `git clang-format` from the repository root. It
  formats the lines that differ from `HEAD` and finds `.clang-format` on its own.
- To format a new file completely, run `clang-format -i <file>`.
- `.clang-format` needs clang-format 14 or newer (it uses `SeparateDefinitionBlocks`).
  Different versions can format some code differently. When the formatter changes lines that
  you did not touch, stop and ask the user which version to use.
- Use `// clang-format off` / `// clang-format on` only for tables that must stay aligned
  (for example variable registration lists).
- Whole-file reformat commits are listed in `.git-blame-ignore-revs`; add the hash there if you
  make one.

### JSON

`exceptions_for_json_formatting.txt` lists the JSON files that are exempt from the JSON
formatting check. The checker itself is not in this repository. Write new JSON files with a
consistent 4-space indentation and do not add them to the exceptions list. When you reformat a
listed file, remove it from the list.

## Check the result

Run `git diff` and make sure that the formatter changed only the files and lines of your
change.
