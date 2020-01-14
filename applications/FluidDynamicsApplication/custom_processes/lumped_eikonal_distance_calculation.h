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

#ifndef KRATOS_LUMPED_EIKONAL_DISTANCE_CALCULATION_H
#define KRATOS_LUMPED_EIKONAL_DISTANCE_CALCULATION_H

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"
#include "custom_utilities/fluid_element_utilities.h"

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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) LumpedEikonalDistanceCalculation : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LumpedEikonalDistanceCalculation
    KRATOS_CLASS_POINTER_DEFINITION(LumpedEikonalDistanceCalculation);

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate paramters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     * @param maxNumIterations Maximum nubmer artificial time-steps to reach the steady-state distance
     * @param tolerance Tolerance in the calculated distance, which determines the steady-state distance
     * @param pseudoTimeStep The size of the artificial time-step to reach the steady-state distance
     */
    LumpedEikonalDistanceCalculation(
        ModelPart& rModelPart,
        const int maxNumIterations,
        const double tolerance,
        const double pseudoTimeStep);

    /// Destructor.
    ~LumpedEikonalDistanceCalculation() override {}

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
     * @brief Initialization of the process 
     */
    void Initialize();

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
        buffer << " LumpedEikonalDistanceCalculation";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << " LumpedEikonalDistanceCalculation";}

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

    // Process parameters
    int mMaxNumIterations = 20;
    double mTolerance = 1.0e-6;
    double mPseudoTimeStep = 1.0;


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This clears the variables
     */
    void ClearVariables();

    /**
     * This gets the auxiliary distance value, which works as distance difference
     * @param rThisGeometry The geometry of the element
     * @param i The node index
     */
    double& GetDDistance(
        Element::GeometryType& rThisGeometry,
        unsigned int i
        );

    /**
     * @brief This divides the distance difference value by the nodal area
     */
    void PonderateDDistance();


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

}; // Class DistanceModificationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_LUMPED_EIKONAL_DISTANCE_CALCULATION_H

