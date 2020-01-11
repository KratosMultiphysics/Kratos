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

#if !defined(KRATOS_FILESYSTEM_AVAILABILITY)
#define KRATOS_FILESYSTEM_AVAILABILITY

// std::filesystem is part of C++17 and not supported by every compiler. Here we check if it is available.
// it is in a separate file due to the order of includes in "includes/kratos_filesystem.h" and "sources/kratos_filesystem.cpp"
#if defined(__cplusplus) && __cplusplus >= 201703L
    #if defined(__has_include) && __has_include(<filesystem>) // has_include is C++17, hence has to be checked in a separatel line
        #define KRATOS_FILESYSTEM_AVAILABLE
    #endif
#endif

#endif // KRATOS_FILESYSTEM_AVAILABILITY defined