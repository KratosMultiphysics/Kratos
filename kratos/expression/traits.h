//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//                   Suneth Warnakulasuriya
//

#pragma once

namespace Kratos {


enum MeshType {
    Local,
    Ghost,
    Interface
}; // enum MeshType

enum ContainerType
{
    NodalHistorical,
    NodalNonHistorical,
    ConditionNonHistorical,
    ElementNonHistorical
}; // enum ContainerType


} // namespace Kratos
