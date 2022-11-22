//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A Roig

#pragma once

// This is a wrapper header to group to handle intel TBB and include
// problems in windows.

// System includes

// External includes
#ifdef KRATOS_INTEL_TBB
    #define NOMINMAX // Fix for windows.h
    #include "tbb/scalable_allocator.h"
#endif

// Project includes