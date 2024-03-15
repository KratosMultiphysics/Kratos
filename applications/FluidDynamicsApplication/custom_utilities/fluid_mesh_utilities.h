//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"

// Application includes


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidMeshUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using ConditionsConnectivityMapType =  std::unordered_map<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> >;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Assigns the conditions' neigbour (parent) elements
     * This function loops the fluid mesh conditions to assign their neighbour (parent) element
     * Note that it also checks for repeated conditions, which are not allowed in fluid meshes
     * @param rModelPart The model part to find the maximum edge length
     * @param CheckRepeatedConditions Checks for repeated conditions
     */
    static void AssignNeighbourElementsToConditions(
        ModelPart& rModelPart,
        const bool CheckRepeatedConditions = true);

    ///@}
};

///@}

} // namespace Kratos

