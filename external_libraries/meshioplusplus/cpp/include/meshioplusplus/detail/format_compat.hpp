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
 * @file format_compat.hpp
 * @brief Portable stand-in for `std::format` on toolchains whose `<format>`
 * is unavailable (e.g. GCC < 13's libstdc++, or clang built against such a
 * libstdc++ - the header does not exist there, so it cannot even be
 * `#include`d, let alone used).
 *
 * Availability is detected through `<version>`'s `__cpp_lib_format` feature
 * test macro rather than `__has_include(<format>)`: `<version>` is
 * guaranteed to exist for any C++20 standard library, and only defines the
 * macro when the library actually implements the feature - so querying it
 * never risks the same "file not found" this header exists to work around.
 *
 * The fallback formatter supports only bare `"{}"` placeholders (no format
 * specs, no positional arguments, no escaping of literal braces) - the only
 * pattern the log/error messages in this codebase use.
 */

// System includes
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <version>

#if !defined(MESHIOPLUSPLUS_FORCE_NO_STD_FORMAT) && defined(__cpp_lib_format) && \
    __cpp_lib_format >= 201907L
#define MESHIOPLUSPLUS_HAS_STD_FORMAT 1
#include <format>
#endif

namespace meshioplusplus {
namespace detail {

#ifdef MESHIOPLUSPLUS_HAS_STD_FORMAT

/** @brief Forwards to `std::format` (compile-time checked format string). */
template <class... Args>
std::string format_compat(std::format_string<Args...> rFmt, Args&&... rArgs) {
    return std::format(rFmt, std::forward<Args>(rArgs)...);
}

#else

/** @brief No-argument overload: the format string, verbatim. */
inline std::string format_compat(std::string_view rFmt) {
    return std::string(rFmt);
}

/**
 * @brief Recursively substitutes each `"{}"` in `rFmt` with the next argument
 * (via `operator<<`), left to right.
 */
template <class T, class... Rest>
std::string format_compat(std::string_view rFmt, const T& rValue, const Rest&... rRest) {
    const std::size_t pos = rFmt.find("{}");
    if (pos == std::string_view::npos)
        return std::string(rFmt);
    std::ostringstream out;
    out << rFmt.substr(0, pos) << rValue;
    return out.str() + format_compat(rFmt.substr(pos + 2), rRest...);
}

#endif

}  // namespace detail
}  // namespace meshioplusplus
