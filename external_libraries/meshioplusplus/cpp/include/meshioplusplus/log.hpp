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
 * @file log.hpp
 * @brief Minimal, header-only logging built on C++20 `std::format` and
 * `std::source_location`.
 *
 * Usage:
 * @code
 * meshioplusplus::log::warn("MED: orientation for '{}' not implemented", type);
 * @endcode
 *
 * Design points:
 *  - Format strings are compile-time checked (`std::format_string`), so a
 *    mismatched `{}` placeholder is a compile error, not a runtime one.
 *  - Every message automatically carries its call site (`file:line`) via a
 *    defaulted `std::source_location` parameter — callers never pass it
 *    explicitly.
 *  - Runtime filtering is controlled by the `MESHIOPLUSPLUS_LOG_LEVEL`
 *    environment variable: `"debug"`, `"info"`, `"warn"` (the default),
 *    `"error"`, or `"off"`. The variable is read exactly once (cached in a
 *    function-local `static`).
 *  - Messages are written to stderr through `std::osyncstream` where the
 *    standard library provides it (guaranteeing concurrent log calls, e.g.
 *    from bodies passed to `parallel_for`, never interleave mid-line); on
 *    standard libraries that ship a `<syncstream>` header without actually
 *    defining `std::osyncstream` (observed with Emscripten's non-threaded
 *    libc++ -- `__cpp_lib_syncbuf` is unset there), a `std::mutex`-guarded
 *    plain write to `std::cerr` gives the same serialization guarantee.
 *  - A filtered-out call costs a single branch: no formatting and no
 *    allocation happen unless the level passes the threshold.
 *  - There is no printf/`std::cerr` logging anywhere else in the codebase;
 *    genuine error conditions remain C++ exceptions (see exceptions.hpp) —
 *    this facility is for diagnostics/warnings only.
 */

// System includes
#include <cstdlib>
#include <format>
#include <iostream>
#include <mutex>
#include <source_location>
#include <string_view>
#include <type_traits>
#include <utility>
#include <version>

#if __has_include(<syncstream>)
#include <syncstream>
#endif

namespace meshioplusplus {
namespace log {

/**
 * @brief Severity levels, ordered so a numerically larger value is louder.
 *
 * `threshold()` returns the minimum level that is actually emitted; a call
 * at a given level is emitted iff `level >= threshold()`. `Off` suppresses
 * every message.
 */
enum class Level : int { Debug = 0, Info = 1, Warn = 2, Error = 3, Off = 4 };

/**
 * @brief The active logging threshold, read once from `MESHIOPLUSPLUS_LOG_LEVEL`.
 *
 * Parses the environment variable on first call (memoized in a function-local
 * `static`, so later changes to the environment have no effect for the
 * lifetime of the process) and defaults to `Level::Warn` when unset or
 * unrecognized. Accepted (case-sensitive) values: `debug`, `info`,
 * `warn`/`warning`, `error`, `off`/`none`.
 *
 * @return The configured minimum `Level` to emit.
 */
inline Level threshold() {
    static const Level lvl = [] {
        const char* env = std::getenv("MESHIOPLUSPLUS_LOG_LEVEL");
        if (env == nullptr)
            return Level::Warn;
        std::string_view s(env);
        if (s == "debug")
            return Level::Debug;
        if (s == "info")
            return Level::Info;
        if (s == "warn" || s == "warning")
            return Level::Warn;
        if (s == "error")
            return Level::Error;
        if (s == "off" || s == "none")
            return Level::Off;
        return Level::Warn;
    }();
    return lvl;
}

/**
 * @brief Whether a message at level `lvl` would actually be emitted.
 * @param lvl The level to test.
 * @return `true` iff `lvl >= threshold()`.
 */
inline bool enabled(Level lvl) {
    return static_cast<int>(lvl) >= static_cast<int>(threshold());
}

/**
 * @brief Formats and writes one log line to stderr.
 *
 * Strips any directory prefix from `loc.file_name()` (keeping just the
 * basename) and writes `"meshio <level> [<file>:<line>] <msg>\n"` through a
 * `std::osyncstream`, which serializes the write against other threads doing
 * the same so concurrent callers (e.g. bodies run under `parallel_for`)
 * never produce interleaved/garbled lines.
 *
 * @param lvl The severity to label the line with.
 * @param msg The already-formatted message body.
 * @param loc The call site to report (normally the caller's, captured via
 *            `FormatWithLocation`).
 */
inline void write(Level lvl, std::string_view msg, const std::source_location& rLoc) {
    constexpr std::string_view names[] = {"debug", "info", "warning", "error"};
    std::string_view file = rLoc.file_name();
    if (auto p = file.find_last_of("/\\"); p != std::string_view::npos)
        file.remove_prefix(p + 1);
    std::string line =
        std::format("meshio {} [{}:{}] {}\n", names[static_cast<int>(lvl)], file, rLoc.line(), msg);
#if defined(__cpp_lib_syncbuf)
    std::osyncstream(std::cerr) << line;
#else
    static std::mutex log_mutex;
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << line;
#endif
}

/**
 * @brief Wraps a compile-time-checked format string with the caller's
 * source location, captured implicitly via a defaulted constructor
 * parameter.
 *
 * This is the trick that lets `log::warn("x={}", x)` both (a) validate the
 * format string against `Args...` at compile time (via `std::format_string`)
 * and (b) automatically know its own call site, without the caller ever
 * writing `std::source_location::current()` themselves: the implicit,
 * `consteval` converting constructor captures `std::source_location::current()`
 * as a default argument evaluated at the *call site* of `debug`/`info`/
 * `warn`/`error`, then packages it alongside the checked format string into
 * one object those functions take by value.
 *
 * @tparam Args The types of the format arguments, used to validate `fmt`.
 */
template <class... Args>
struct FormatWithLocation {
    std::format_string<Args...> mFmt;
    std::source_location mLoc;

    template <class S>
    consteval FormatWithLocation(  // NOLINT(google-explicit-constructor)
        const S& rS, std::source_location l = std::source_location::current())
        : mFmt(rS), mLoc(l) {}
};

/**
 * @brief Logs a debug-level message (lowest severity; off by default).
 *
 * No-op (no formatting, no allocation) unless `MESHIOPLUSPLUS_LOG_LEVEL=debug`.
 * @tparam Args Format argument types, deduced from `args`.
 * @param f Compile-time-checked format string + implicit call site.
 * @param args Values substituted into `f`.
 */
template <class... Args>
void debug(FormatWithLocation<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Debug))
        return;
    write(Level::Debug, std::format(f.mFmt, std::forward<Args>(args)...), f.mLoc);
}

/**
 * @brief Logs an info-level message.
 *
 * No-op unless `MESHIOPLUSPLUS_LOG_LEVEL` is `debug` or `info`.
 * @tparam Args Format argument types, deduced from `args`.
 * @param f Compile-time-checked format string + implicit call site.
 * @param args Values substituted into `f`.
 */
template <class... Args>
void info(FormatWithLocation<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Info))
        return;
    write(Level::Info, std::format(f.mFmt, std::forward<Args>(args)...), f.mLoc);
}

/**
 * @brief Logs a warn-level message. This is the default active threshold.
 *
 * Used, for example, when a C++ format implementation encounters a
 * recognized-but-unhandled construct and degrades gracefully rather than
 * failing (e.g. an unimplemented MED orientation).
 * @tparam Args Format argument types, deduced from `args`.
 * @param f Compile-time-checked format string + implicit call site.
 * @param args Values substituted into `f`.
 */
template <class... Args>
void warn(FormatWithLocation<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Warn))
        return;
    write(Level::Warn, std::format(f.mFmt, std::forward<Args>(args)...), f.mLoc);
}

/**
 * @brief Logs an error-level message.
 *
 * @note This is diagnostic logging only, not the mechanism for signaling
 * failures to callers — genuine error conditions must still be reported by
 * throwing `ReadError`/`WriteError` (see exceptions.hpp); this call alone
 * does not stop execution or propagate anything.
 * @tparam Args Format argument types, deduced from `args`.
 * @param f Compile-time-checked format string + implicit call site.
 * @param args Values substituted into `f`.
 */
template <class... Args>
void error(FormatWithLocation<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Error))
        return;
    write(Level::Error, std::format(f.mFmt, std::forward<Args>(args)...), f.mLoc);
}

}  // namespace log
}  // namespace meshioplusplus
