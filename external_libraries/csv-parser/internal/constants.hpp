/** @file
 *  Defines CSV global constants
 */

#pragma once
#include <deque>

#include "csv_format.hpp"

#if defined(_WIN32)
#include <Windows.h>
#undef max
#undef min
#elif defined(__linux__)
#include <unistd.h>
#endif

namespace csv {
    namespace internals {
        // PAGE_SIZE macro could be already defined by the host system.
        #if defined(PAGE_SIZE)
        #undef PAGE_SIZE
        #endif

        // Get operating system specific details
        #if defined(_WIN32)
            inline int getpagesize() {
                _SYSTEM_INFO sys_info = {};
                GetSystemInfo(&sys_info);
                return sys_info.dwPageSize;
            }

            /** Size of a memory page in bytes */
            const int PAGE_SIZE = getpagesize();
        #elif defined(__linux__) 
            const int PAGE_SIZE = getpagesize();
        #else
            const int PAGE_SIZE = 4096;
        #endif

        /** For functions that lazy load a large CSV, this determines how
         *  many bytes are read at a time
         */
        constexpr size_t ITERATION_CHUNK_SIZE = 50000000; // 50MB
    }

    /** Integer indicating a requested column wasn't found. */
    constexpr int CSV_NOT_FOUND = -1;

    /** Used for counting number of rows */
    using RowCount = long long int;
}