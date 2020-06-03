/** @file
 *  Defines various compatibility macros
 */

/** Helper macro which should be #defined as "inline" 
 *  in the single header version
 */
#define CSV_INLINE

#pragma once
#include "../external/string_view.hpp"

// If there is another version of Hedley, then the newer one 
// takes precedence.
// See: https://github.com/nemequ/hedley
#include "../external/hedley.h"

namespace csv {
    /**
     *  @def IF_CONSTEXPR
     *  Expands to `if constexpr` in C++17 and `if` otherwise
     *
     *  @def CONSTEXPR_VALUE
     *  Expands to `constexpr` in C++17 and `const` otherwise.
     *  Mainly used for global variables.
     *
     *  @def CONSTEXPR
     *  Expands to `constexpr` in C++17 and `inline` otherwise.
     *  Intended for functions and methods.
     */

    #if CMAKE_CXX_STANDARD == 17 || __cplusplus >= 201703L
        #define CSV_HAS_CXX17
    #endif

    #ifdef CSV_HAS_CXX17
        #include <string_view>
        /** @typedef string_view
         *  The string_view class used by this library.
         */
        using string_view = std::string_view;
    #else
        /** @typedef string_view
         *  The string_view class used by this library.
         */
        using string_view = nonstd::string_view;
    #endif

    #ifdef CSV_HAS_CXX17
        #define IF_CONSTEXPR if constexpr
        #define CONSTEXPR_VALUE constexpr
    #else
        #define IF_CONSTEXPR if
        #define CONSTEXPR_VALUE const
    #endif

    // Resolves g++ bug with regard to constexpr methods
    #if defined __GNUC__ && !defined __clang__
        #if __GNUC__ >= 7
            #if defined(CSV_HAS_CXX17) && (__GNUC_MINOR__ >= 2 || __GNUC__ >= 8)
                #define CONSTEXPR constexpr
            #endif
        #endif
    #else
        #ifdef CSV_HAS_CXX17
            #define CONSTEXPR constexpr
        #endif
    #endif

    #ifndef CONSTEXPR
        #define CONSTEXPR inline
    #endif
}
