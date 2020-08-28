//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    ME
//
//

#ifndef KRATOS_HYPERBOLIC_DISTANCE_H
#define KRATOS_HYPERBOLIC_DISTANCE_H

// Just to follow the Semi-Lagrangian code templates
#define PRESSURE_ON_EULERIAN_MESH
#define USE_FEW_PARTICLES

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "utilities/binbased_fast_point_locator.h"
#include "processes/compute_nodal_gradient_process.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Utility to due the redistancing based on the time-dependednt Eikonal equation

    constexpr double mesh_tolerance = 1.0e-9;

template<std::size_t TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) HyperbolicDistanceReinitialization : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HyperbolicDistanceReinitialization
    KRATOS_CLASS_POINTER_DEFINITION(HyperbolicDistanceReinitialization<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate paramters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     * @param pSearchStructure Pointer to utility for searching a specific position in mesh
     * @param pComputeGradient Pointer to process for calculating nodal gradient
     * @param maxNumIterations Maximum nubmer artificial time-steps to reach the steady-state distance
     * @param tolerance Tolerance in the calculated distance, which determines the steady-state distance
     * @param pseudoTimeStep The size of the artificial time-step to reach the steady-state distance
     */
    HyperbolicDistanceReinitialization(
        ModelPart& rModelPart,
        typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure,
        ComputeNodalGradientProcess<true>::Pointer pComputeGradient,
        const int maxNumIterations = 20,
        const double tolerance = 1.0e-6,
        const double pseudoTimeStep = 1.0);

    /// Destructor.
    ~HyperbolicDistanceReinitialization() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Execution of the process
     */
    void Execute();

    // ///@}
    // ///@name Inquiry
    // ///@{

    // ///@}
    // ///@name Input and output
    // ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << " HyperbolicDistanceReinitialization";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << " HyperbolicDistanceReinitialization";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // Reference to the model part
    const ModelPart& mrModelPart;
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;
    ComputeNodalGradientProcess<true>::Pointer mpComputeGradient;

    // Process parameters
    int mMaxNumIterations;
    double mTolerance;
    double mPseudoTimeStep;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class HyperbolicDistanceReinitialization

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_HYPERBOLIC_DISTANCE_H

