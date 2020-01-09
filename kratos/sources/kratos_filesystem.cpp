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

// Project includes
#include "utilities/check_filesystem_availability.h" // has to be included before using "KRATOS_FILESYSTEM_AVAILABLE"

// External includes
#ifndef KRATOS_FILESYSTEM_AVAILABLE
// only required if we ghc/filesystem is used
// has to be included before "includes/kratos_filesystem.h", see following explanation

//---------------------------------------------------------------------------------------
// fs_impl.hpp - The implementation header for the header/implementation seperated usage of
//               ghc::filesystem.
// This file can be used to hide the implementation of ghc::filesystem into a single cpp.
// The cpp has to include this before including fs_fwd.hpp directly or via a different
// header to work.
//---------------------------------------------------------------------------------------
#include "ghc/fs_impl.hpp"
#endif

// Project includes
#include "includes/kratos_filesystem.h"
