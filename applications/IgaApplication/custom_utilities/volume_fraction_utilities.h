//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) VolumeFractionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the active volume fraction of every knot span of an axis-aligned 2D background patch.
     *
     * @details The domain is the region inside the outer skin loops (or the whole span box if the outer skin has
     * no conditions) minus the region inside the inner skin loops. Loops are chained by node connectivity and may
     * have any orientation.
     *
     * @param rOuterSkinModelPart Closed loops of 2-noded line conditions bounding the domain from outside. May be empty.
     * @param rInnerSkinModelPart Closed loops of 2-noded line conditions bounding the holes of the domain. May be empty.
     * @param rSpansX Strictly increasing span boundaries in x, in the same coordinates as the skin nodes.
     * @param rSpansY Strictly increasing span boundaries in y, in the same coordinates as the skin nodes.
     * @return Matrix with one entry per knot span: entry (i, j) is the volume fraction of the i-th span in x
     * and the j-th span in y.
     */
    static Matrix ComputeKnotSpanVolumeFractions(
        const ModelPart& rOuterSkinModelPart,
        const ModelPart& rInnerSkinModelPart,
        const std::vector<double>& rSpansX,
        const std::vector<double>& rSpansY);

    ///@}
};

///@}

} // namespace Kratos
