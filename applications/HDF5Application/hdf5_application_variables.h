//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

// Application includes
#include "custom_utilities/mesh_location_container.h"

namespace Kratos
{

    KRATOS_DEFINE_APPLICATION_VARIABLE(HDF5_APPLICATION, HDF5::MeshLocationContainer::Pointer, HDF5_MESH_LOCATION_CONTAINER)

}  // namespace Kratos.