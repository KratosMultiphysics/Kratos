//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//

#ifndef KRATOS_DERIVATIVE_RECOVERY_PROCESS_H
#define KRATOS_DERIVATIVE_RECOVERY_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{
///@addtogroup SwimmingDEMApplication
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

class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DerivativeRecoveryProcess
    KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecoveryProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    DerivativeRecoveryProcess(
        ModelPart& rModelPart,
        Parameters Parameters);

    /// Constructor with Kratos model
    DerivativeRecoveryProcess(
        Model& rModel,
        Parameters Parameters);

    /// Destructor.
    ~DerivativeRecoveryProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Operations
    ///@{

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
        buffer << "DerivativeRecoveryProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "DerivativeRecoveryProcess";}

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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters rParameters);

    void CalculateScalarGradients();

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
    DerivativeRecoveryProcess() = delete;

    /// Assignment operator.
    DerivativeRecoveryProcess& operator=(DerivativeRecoveryProcess const& rOther) = delete;

    /// Copy constructor.
    DerivativeRecoveryProcess(DerivativeRecoveryProcess const& rOther) = delete;

    ///@}

}; // Class DerivativeRecoveryProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DERIVATIVE_RECOVERY_PROCESS_H
