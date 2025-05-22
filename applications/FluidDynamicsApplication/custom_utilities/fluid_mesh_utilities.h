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
     * @brief Checks if elements are simplex
     * This function loops all the elements in the provided model part to check if all of
     * them are simplex (simplicial complex), that is, linear triangles and tetrahedra
     * @param rModelPart The model part to which the check is performed
     * @return true All elements are siplex
     * @return false There are non-simplex elements
     */
    static bool AllElementsAreSimplex(const ModelPart& rModelPart);

    /**
     * @brief Assigns the conditions' neighbour (parent) elements
     * This function loops the fluid mesh conditions to assign their neighbour (parent) element
     * Note that it also checks for repeated conditions, which are not allowed in fluid meshes
     * Also note that this function assumes that elements' face geometries match the conditions' ones
     * @param rModelPart The model part to which the parents are to be assigned
     * @param CheckRepeatedConditions Checks for repeated conditions
     */
    static void AssignNeighbourElementsToConditions(
        ModelPart& rModelPart,
        const bool CheckRepeatedConditions = true);

    ///@}
};

///@}

} // namespace Kratos

