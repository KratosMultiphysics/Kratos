//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/container_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ImplicitFilterUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static void CalculateNodeNeighbourCount(
        ModelPart& rModelPart);

    static void SetBulkRadiusForShapeFiltering(
        ModelPart& rModelPart);

    static void AssignVectorNodalExpressionToScalarVariable(const Variable<double>& rVariable,
        const SpecializedContainerExpression<ModelPart::NodesContainerType,
        ContainerDataIO<ContainerDataIOTags::NonHistorical>>& rContainer, int step);

    static void AssignScalarVariableToVectorNodalExpression(const Variable<double>& rVariable,
        SpecializedContainerExpression<ModelPart::NodesContainerType,
        ContainerDataIO<ContainerDataIOTags::NonHistorical>>& rContainer, int step);

    ///@}
};

///@}
}