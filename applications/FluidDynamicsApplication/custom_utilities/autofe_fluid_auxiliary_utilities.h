//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
// #include "modified_shape_functions/modified_shape_functions.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) AutoFEFluidAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Postprocess the midpoint nodes for a given scalar variable in P2P1 elements
     * This function takes the edges' midpoint nodes in P2P1 elements and postprocess the given variable, which
     * is assumed to be stored as a historical variable, from the edges' endpoint values.
     * Note that the nodal flag VISITED is used to mark the nodes which pressure has been already set.
     * @param rModelPart The model part to which the variable is to be postprocessed
     * @param rVariable The scalar variable to be processed
     */
    static void PostprocessP2P1ContinuousScalarVariable(ModelPart& rModelPart, const Variable<double>& rVariable);

    ///@}
private:

    /**
     * @brief Auxilary function to postprocess one P2P1 edge variable
     * This function postprocesses a given scalar valriable in a P2P1 element edge midpoint.
     * Once the variable's value is set, the node is marked as VISITED.
     * @param rGeometry Reference to current element geometry
     * @param rVarialbe Reference to the variable to process
     * @param PostNodeLocalId Local id of the node to which the varaible value is to be set
     * @param EdgeNodeLocalIdI Local id of the i-node of the edge to which previous node belongs
     * @param EdgeNodeLocalIdJ Local id of the j-node of the edge to which previous node belongs
     */
    static void PostprocessP2P1NodeScalarVariable(
        GeometryType& rGeometry,
        const Variable<double>& rVariable,
        const std::size_t PostNodeLocalId,
        const std::size_t EdgeNodeLocalIdI,
        const std::size_t EdgeNodeLocalIdJ)
    {
        if (rGeometry[PostNodeLocalId].IsNot(VISITED)) {
            rGeometry[PostNodeLocalId].SetLock();
            const double p_i = rGeometry[EdgeNodeLocalIdI].FastGetSolutionStepValue(rVariable);
            const double p_j = rGeometry[EdgeNodeLocalIdJ].FastGetSolutionStepValue(rVariable);
            rGeometry[PostNodeLocalId].FastGetSolutionStepValue(rVariable) = 0.5 * (p_i + p_j);
            rGeometry[PostNodeLocalId].Set(VISITED, true);
            rGeometry[PostNodeLocalId].UnSetLock();
        }
    }

};

///@}

} // namespace Kratos
