//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <map>
#include "geometries/geometry_data.h"

namespace Kratos {

class VtkDefinitions {
public:
    ///@name Public static member variables
    ///@{

    // IMPORTANT: The map KratosVtkGeometryTypes is to be extended to support new geometries
    // NOTE: See https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf and https://vtk.org/wp-content/uploads/2021/08/VTKUsersGuide.pdf
    static const std::map<GeometryData::KratosGeometryType, char> KratosVtkGeometryTypes;

    ///@}
};

} // namespace Kratos