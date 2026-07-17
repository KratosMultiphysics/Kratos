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
 * @file source_location_compat.hpp
 * @brief Portable stand-in for `std::source_location` on toolchains whose
 * `<source_location>` does not actually populate `std::source_location`
 * (observed with clang-14 against some libstdc++ versions: the header
 * includes without error, but `std::source_location` is simply not declared).
 *
 * Availability is detected through `<version>`'s `__cpp_lib_source_location`
 * feature test macro. The fallback is implemented with the same
 * `__builtin_FILE()`/`__builtin_LINE()` compiler builtins the standard
 * implementations themselves are built on (supported by both GCC and Clang),
 * so captured call sites are identical to the real thing.
 */

// System includes
#include <cstdint>
#include <version>

#if defined(__cpp_lib_source_location) && __cpp_lib_source_location >= 201907L
#define MESHIOPLUSPLUS_HAS_STD_SOURCE_LOCATION 1
#include <source_location>
#endif

namespace meshioplusplus {
namespace detail {

#ifdef MESHIOPLUSPLUS_HAS_STD_SOURCE_LOCATION

using source_location = std::source_location;

#else

class source_location {
public:
    static consteval source_location current(
        const char* pFile = __builtin_FILE(), int Line = __builtin_LINE()) noexcept {
        source_location loc;
        loc.mFile = pFile;
        loc.mLine = Line;
        return loc;
    }

    constexpr const char* file_name() const noexcept { return mFile; }
    constexpr std::uint_least32_t line() const noexcept { return static_cast<std::uint_least32_t>(mLine); }

private:
    const char* mFile = "";
    int mLine = 0;
};

#endif

}  // namespace detail
}  // namespace meshioplusplus
