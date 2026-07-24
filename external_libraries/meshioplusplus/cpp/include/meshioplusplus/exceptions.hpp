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
 * @file exceptions.hpp
 * @brief meshio I/O exception types thrown by the C++ core's readers/writers.
 *
 * These are the only exception types the C++ format readers/writers throw on
 * I/O failure (malformed input, unsupported constructs, filesystem errors,
 * etc.). The pybind11 binding layer catches them and re-raises the
 * equivalent Python `meshioplusplus.ReadError` / `meshioplusplus.WriteError`
 * classes, so callers on the Python side see identical behaviour whether a
 * format is handled by the C++ core or by the pure-Python fallback. Because
 * the shim pattern (`__init__.py`) catches *any* exception from the C++ path
 * to decide whether to fall back to Python, throwing these (rather than
 * e.g. asserting or returning error codes) is what makes that fallback work.
 */

// System includes
#include <stdexcept>
#include <string>

namespace meshioplusplus {

/**
 * @brief Thrown by C++ readers when the input file/stream cannot be parsed.
 *
 * Covers malformed content, missing required sections, and unsupported
 * constructs that a given format's C++ reader deliberately does not handle
 * (in which case the format's Python shim catches this and falls back to the
 * pure-Python reference reader). Maps 1:1 to Python's `meshioplusplus.ReadError`.
 */
struct ReadError : std::runtime_error {
    ReadError() : std::runtime_error("") {}
    explicit ReadError(const std::string& rMsg) : std::runtime_error(rMsg) {}
};

/**
 * @brief Thrown by C++ writers when a mesh cannot be serialized to a format.
 *
 * Covers unsupported cell types, ragged/ill-formed mesh data the writer does
 * not accept, and any other output-side constraint violation (in which case
 * the format's Python shim catches this and falls back to the pure-Python
 * reference writer). Maps 1:1 to Python's `meshioplusplus.WriteError`.
 */
struct WriteError : std::runtime_error {
    WriteError() : std::runtime_error("") {}
    explicit WriteError(const std::string& rMsg) : std::runtime_error(rMsg) {}
};

}  // namespace meshioplusplus
