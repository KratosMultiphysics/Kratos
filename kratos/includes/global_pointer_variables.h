//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_components.h"
#include "containers/global_pointers_vector.h"
#include "includes/node.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

namespace Kratos
{
    KRATOS_DEFINE_VARIABLE(GlobalPointersVector<Node>, NEIGHBOUR_NODES)
    KRATOS_DEFINE_VARIABLE(GlobalPointersVector<Node>, NEIGHBOUR_CONDITION_NODES)
    KRATOS_DEFINE_VARIABLE(GlobalPointersVector<Node>, FATHER_NODES)
}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT