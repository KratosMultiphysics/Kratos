//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Bastian Devresse,
//

#ifndef SUBDIVISION_UTILITIES_H
#define SUBDIVISION_UTILITIES_H

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

class KRATOS_API(OPTIMIZATION_APPLICATION) SDSUtils
{

public:
    ///@name Type definitions
    ///@{

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{


    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectForward(
        const ContainerExpression<TContainerType>& rInputExpression);

    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectBackward(
        const ContainerExpression<TContainerType>& rInputExpression);

    static std::vector<double> CalculateMappingRelation(
        ModelPart& rControlPolygon,
        const ModelPart& rControlledMesh,
        const std::string SubdivScheme,
        const bool FixFreeEdges);

    ///@}
};

///@}
}

#endif // SUBDIVISION_UTILITIES_H defined
