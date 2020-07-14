//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

#ifndef KRATOS_ONE_SIDE_PRESSURE_GRADIENT_H
#define KRATOS_ONE_SIDE_PRESSURE_GRADIENT_H

// Project includes
#include "processes/process.h"


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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputeOneSideNodalPressureGradientProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LumpedInterfacePositiveNegativePressureGradient
    KRATOS_CLASS_POINTER_DEFINITION(ComputeOneSideNodalPressureGradientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate paramters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     */
    ComputeOneSideNodalPressureGradientProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    ComputeOneSideNodalPressureGradientProcess(
        ModelPart& rModelPart,
        Parameters Parameters);

    /// Constructor with Kratos model
    ComputeOneSideNodalPressureGradientProcess(
        Model& rModel,
        Parameters Parameters);

    /// Destructor.
    ~ComputeOneSideNodalPressureGradientProcess() override {}

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
    void Execute() override;

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
        buffer << "ComputeOneSideNodalPressureGradientProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "ComputeOneSideNodalPressureGradientProcess";}

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
    ModelPart& mrModelPart;

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

}; // Class ComputeOneSideNodalPressureGradientProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_ONE_SIDE_PRESSURE_GRADIENT_H


