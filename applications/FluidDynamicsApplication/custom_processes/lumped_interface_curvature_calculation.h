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

#ifndef KRATOS_LUMPED_INTERFACE_CURVATURE_CALCULATION_H
#define KRATOS_LUMPED_INTERFACE_CURVATURE_CALCULATION_H

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
//#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) LumpedInterfaceCurvatureCalculation : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LumpedInterfaceCurvatureCalculation
    KRATOS_CLASS_POINTER_DEFINITION(LumpedInterfaceCurvatureCalculation);

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate paramters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     */
    LumpedInterfaceCurvatureCalculation(
        ModelPart& rModelPart);

    /// Destructor.
    ~LumpedInterfaceCurvatureCalculation() override {}

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
        buffer << " LumpedInterfaceCurvatureCalculation";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << " LumpedInterfaceCurvatureCalculation";}

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

}; // Class LumpedInterfaceCurvatureCalculation

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_LUMPED_INTERFACE_CURVATURE_CALCULATION_H


