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
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_HYDRAULIC_FLUID_AUXILIARY_UTILITIES_H)
#define KRATOS_HYDRAULIC_FLUID_AUXILIARY_UTILITIES_H

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "includes/global_pointer_variables.h"

// Application includes
#include "hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) HydraulicFluidAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    ///@}
    ///@name Static Operations
    ///@{
    /**
     * @brief This functions calculates the wetted perimeter of a given condition in order to apply the
     * corresponding inlet boundary condition.
     * @param rModelPart Fluid Model Part
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return 
     */
    static double CalculateWettedPetimeter(
        ModelPart &rModelPart,
        const Flags &rSkinFlag);
    ///@}
    ///@{
    /**
     * @brief This functions calculates the inlet discharge assuming a critical depth, Froude Number of 1. 
     * @param rModelPart The model part to find the maximum edge length
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return
     */
    static double CalculateInletDischarge(
        ModelPart &rModelPart,
        const Flags &rSkinFlag);
    ///@}
    ///@{
    /**
     * @brief This functions calculates the inlet discharge assuming a critical depth, Froude Number of 1.
     * @param rModelPart The model part to find the maximum edge length
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return
     */
    static double CalculateWettedArea(
        ModelPart &rModelPart,
        const Flags &rSkinFlag);
    ///@}

    ///@{
    /**
     * @brief This functions calculates the inlet discharge assuming a critical depth, Froude Number of 1.
     * @param rModelPart The model part to find the maximum edge length
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return
     */
    static double CalculateConditionArea(
        const GeometryType &rGeometry);
    ///@}

private:
    
    ///@}

}; // namespace Kratos
}
#endif // KRATOS_HYDRAULIC_FLUID_AUXILIARY_UTILITIES_H
