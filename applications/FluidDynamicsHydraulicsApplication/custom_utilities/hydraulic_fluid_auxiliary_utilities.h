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

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/global_pointer_variables.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "../../FluidDynamicsApplication/custom_utilities/fluid_auxiliary_utilities.h"

// Application includes
#include "hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsHydraulicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_HYDRAULICS_APPLICATION) HydraulicFluidAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PointsArrayType = typename GeometryType::PointsArrayType;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    ///@}
    ///@name Static Operations
    ///@{
    /**
     * @brief This functions calculates the wetted perimeter of a given condition in order to apply the
     * corresponding inlet boundary condition.
     * @param rModelPart Inlet Model Part
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param rDistanceVariable Reference to the variable containing the distance
     * @param IsHistorical True if the distance is in the historical database, false otherwise
     * @return Wetted perimeter 
     */
    static double CalculateWettedPetimeter(
        ModelPart &rModelPart,
        const Flags &rSkinFlag,
        const Variable<double> &rDistanceVariable,
        const bool IsHistorical);
    ///@}
    ///@{
    /**
     * @brief This functions calculates the wetted area of a given condition in order to apply the
     * corresponding inlet boundary condition.
     * @param rModelPart Inlet Model Part
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param rDistanceVariable Reference to the variable containing the distance
     * @param IsHistorical True if the distance is in the historical database, false otherwise
     * @return Wetted area
     */
    static double CalculateWettedArea(
        ModelPart &rModelPart,
        const Flags &rSkinFlag,
        const Variable<double> &rDistanceVariable,
        const bool IsHistorical);

    /**
     * @brief Calculates initial water depth guess by taking the average between the maximum and minimum coordinates.
     *@param rModelPart Inlet Model Part
     * @return Initial water depth guess
     */
    static double InitialWaterDepth(ModelPart &rModelPart);

    /**
     * @brief
     * Assign the inlet velocity to all nodes that are wet in the input model part. For dry nodes(air)is assumed that the inlet velocity is null.
     * @param  rModelPart Inlet Model Part
     * @param  InletVelocity the velocity value to be assigned to wet nodes.
     * @param  rDistancesVariable Variable name of the inlet distance.
     */
    static void SetInletVelocity(ModelPart &rModelPart, double InletVelocity, const Variable<double> &rDistanceVariable);

    /**
     * @brief Free the inlet velocity in the nodes belonging to inlet model part.
     * @param  rModelPart Inlet Model Part
     */
    static void FreeInlet(ModelPart& rModelPart);

    /**
     * @brief  Set the free surface (DISTANCE) in the rModelPart equal to the water depth corresponding to Froude 1
     * @param  rModelPart Inlet Model Part
     * @param  rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param  rDistancesVariable Variable name of the inlet distance.
     */
    static void SetInletFreeSurface(ModelPart &rModelPart, const Flags &rSkinFlag, const Variable<double> &rDistanceVariable);

    /**
     * @brief  Artificial viscosity is calculated. The purpose of adding this artificial viscosity is to avoid non-physical spikes in velocities.
     * @param  rModelPart Fluid Model Part
     * @param  WaterDynamicViscosityMax It is a threshold value to prevent adding excessive artificial numerical viscosity and thereby losing the real physics.
     */
    static void CalculateArtificialViscosity(ModelPart &rModelPart, double WaterDynamicViscosityMax);

    /**
     * @brief  When there is inflow on a boundary considered as an outlet, this function retains only the tangential component, preventing inflows that cause instabilities.
     * @param  rModelPart Fluid Model Part
     * @param  rVariable it possible to use the variable VELOCITY_FRACTIONAL or VELOCITY
     */

    static void ApplyOutletInflowLimiter(ModelPart &rModelPart,const Variable<array_1d<double, 3>>& rVariable);

     

    ///@}

private :

    struct EdgeDataContainer
    {
        NodeType::Pointer pNodeI = nullptr;
        NodeType::Pointer pNodeJ = nullptr;
        SizeType NumberOfRepetitions = 0;
    };

}; // namespace Kratos
}
