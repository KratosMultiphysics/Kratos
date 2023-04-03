//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    clabra
//

#pragma once

// External includes

#if defined(KRATOS_INDEPENDENT)

// clean definitions of kratos macros
#define KRATOS_CLASS_POINTER_DEFINITION(variable) \
  typedef variable* Pointer

#define KRATOS_WATCH(variable)

#else

// include kratos definitions
#include "includes/define.h"

#endif // KRATOS_INDEPENDENT_LIBRARY

// Project includes
#include "tree.h"
#include "bucket.h"
#include "kd_tree.h"
#include "octree.h"
#include "octree_binary.h"
#include "bins_static.h"
#include "bins_dynamic.h"
#include "bins_dynamic_objects.h"
#include "bins_static_objects.h"