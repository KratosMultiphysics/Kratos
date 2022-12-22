//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
//

#if !defined(KRATOS_SHOCK_CAPTURING_UTILITIES_H)
#define KRATOS_SHOCK_CAPTURING_UTILITIES_H

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "modified_shape_functions/modified_shape_functions.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template<std::size_t TDim, std::size_t TNumNodes>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Static Operations
    ///@{

    static void ExecuteRMethodType(const ModelPart &mrModelPart);

    /**
     * @brief Calculate the fluid negative volume
     * For the given model part, this function calculates the total negative fluid volume fraction.
     * It is assumed that a unique element geometry type it is present in the mesh.
     * Only simplex geometries (linear triangle and tetrahedron) are supported.
     * The negative fraction must be described in terms of a continuous level set function stored
     * in the variable DISTANCE of the historical database
     * @param rModelPart The model part to calculate the negative volume
     * @return double Fluid negative volume
     */

};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_AUXILIARY_UTILITIES_H
