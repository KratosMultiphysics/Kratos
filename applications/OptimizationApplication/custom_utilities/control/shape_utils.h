//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ShapeUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static void GetNodalCoordinates(
        ContainerExpression<ModelPart::NodesContainerType>& rInput);

    static void SetNodalCoordinates(
        ContainerExpression<ModelPart::NodesContainerType>& rInput);

    static void UpdateNodalCoordinates(
        ContainerExpression<ModelPart::NodesContainerType>& rInput);

    ///@}
};

///@}
}