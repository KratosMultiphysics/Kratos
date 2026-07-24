//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file registry.hpp
 * @brief C++-level format dispatch registry shared by the flat bindings
 *        (WASM/JS in `bindings_js/`, the C API in `bindings_c/`).
 *
 * The pybind11 binding (`bindings/_core.cpp`) exposes one function per format
 * and leaves extension dispatch entirely to Python (`_helpers.py`); the flat
 * bindings instead need a C++-side `format name -> read/write function` table
 * plus an `extension -> default format` map. Those tables originally lived in
 * `bindings_js/js_bindings.cpp`; they are hoisted here (compiled into
 * `meshioplusplus_core_obj` via `cpp/src/registry.cpp`) so the JS and C
 * bindings share one copy that cannot drift.
 *
 * Parameterized writers get a fixed default here (documented per entry in
 * registry.cpp, matching each format's own Python reference default);
 * per-call overrides are a possible future API addition, deliberately out of
 * scope for v1 of both flat bindings.
 *
 * HDF5-backed formats (cgns, h5m, hmf, med, plus XDMF's HDF data path) and
 * netCDF-backed ones (exodus) are registered only when the corresponding
 * `MESHIOPLUSPLUS_HAS_*` macro is defined -- never under Emscripten, so the
 * WASM format set is unchanged by this refactor. Their extensions are mapped
 * unconditionally so that a build without the dependency reports "format
 * compiled out" (see registry_compiled_out()) instead of the misleading
 * "cannot infer format".
 */

// System includes
#include <functional>
#include <map>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

using ReadFn = std::function<Mesh(const std::string&)>;
using WriteFn = std::function<void(const std::string&, const Mesh&)>;

/** @brief `format name -> reader` for every format readable in this build. */
const std::map<std::string, ReadFn>& registry_readers();

/** @brief `format name -> writer` for every format writable in this build
 *         (read-only formats like `openfoam` have no entry). */
const std::map<std::string, WriteFn>& registry_writers();

/**
 * @brief `extension (with leading dot) -> default format name`.
 *
 * Ambiguous extensions get this repo's own import-order default (`.msh` ->
 * gmsh, `.inp` -> abaqus); pass an explicit format to select ansys/freefem
 * (.msh) or ansysinp (.inp) instead. Extensions of optional-dependency
 * formats are present even when the format itself is compiled out.
 */
const std::map<std::string, std::string>& registry_extension_defaults();

/**
 * @brief Resolve the effective format: `rFormat` if non-empty, else the
 *        extension default for `rPath`.
 * @throws ReadError if `rFormat` is empty and the extension is unknown.
 */
std::string resolve_format(const std::string& rPath, const std::string& rFormat);

/**
 * @brief The optional dependency a known-but-absent format was compiled out
 *        with, or `nullptr`.
 * @return `"HDF5"` / `"netCDF"` when `rFormat` names a format this build
 *         excluded for lack of that dependency; `nullptr` for formats that
 *         are present or simply unknown. Lets bindings say "format 'med' is
 *         not available in this build (requires HDF5)" instead of "unknown
 *         format".
 */
const char* registry_compiled_out(const std::string& rFormat);

}  // namespace meshioplusplus
