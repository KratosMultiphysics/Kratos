//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_FILESYSTEM)
#define  KRATOS_FILESYSTEM

// System / External includes
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__has_include) && __has_include(<filesystem>)
    #include <filesystem>
    #define KRATOS_FILESYSTEM_AVAILABLE
#else
    #include "ghc/fs_fwd.hpp"
#endif

namespace Kratos {

// alias to make it available as Kratos::filesystem
#ifdef KRATOS_FILESYSTEM_AVAILABLE
    namespace filesystem = std::filesystem;
#else
    namespace filesystem = ghc::filesystem;
#endif
}

#endif // KRATOS_FILESYSTEM defined