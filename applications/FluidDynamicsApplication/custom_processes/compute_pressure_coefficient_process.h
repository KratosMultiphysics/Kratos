//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

#ifndef KRATOS_COMUTE_PRESSURE_COEFFICIENT_PROCESS_H
#define KRATOS_COMUTE_PRESSURE_COEFFICIENT_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

// Application includes


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

/// This process computes the pressure coefficient as a function of reference fluid properties
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputePressureCoefficientProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    using NodeType = ModelPart::NodeType;

    /// Pointer definition of ComputePressureCoefficientProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputePressureCoefficientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos model
    ComputePressureCoefficientProcess(
        Model& rModel,
        Parameters Params);

    /// Destructor.
    ~ComputePressureCoefficientProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    const Parameters GetDefaultParameters() const override;

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    void ExecuteBeforeOutputStep() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ComputePressureCoefficientProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "ComputePressureCoefficientProcess";}

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

    ModelPart& mrModelPart;
    bool mComputeAsPostProcess;
    double mFreestreamStaticPressure;   // Freestream pressure
    double mFreestreamDynamicPressure;  // Freestream q=rho*V²/2
    std::function<double(const NodeType&)> mGetPressure;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    void SelectExecutionTime(Parameters Params);

    void SelectPressureGetter(Parameters Params);

    void ReadFreestreamValues(Parameters Params);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    ComputePressureCoefficientProcess() = delete;

    /// Assignment operator.
    ComputePressureCoefficientProcess& operator=(ComputePressureCoefficientProcess const& rOther) = delete;

    /// Copy constructor.
    ComputePressureCoefficientProcess(ComputePressureCoefficientProcess const& rOther) = delete;

    ///@}

}; // Class ComputePressureCoefficientProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_COMUTE_PRESSURE_COEFFICIENT_PROCESS_H
