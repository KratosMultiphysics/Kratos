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

#ifndef KRATOS_STANDARD_RECOVERY_PROCESS_H
#define KRATOS_STANDARD_RECOVERY_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "derivative_recovery_process.h"
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

class KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryProcess : public DerivativeRecoveryProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StandardRecoveryProcess
    KRATOS_CLASS_POINTER_DEFINITION(StandardRecoveryProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor with Kratos parameters.
    StandardRecoveryProcess(
        ModelPart& rModelPart,
        Parameters Param);

    /// Constructor with Kratos model
    StandardRecoveryProcess(
        Model& rModel,
        Parameters Param);

    /// Destructor.
    ~StandardRecoveryProcess() override {}

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
        buffer << "StandardRecoveryProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "StandardRecoveryProcess";}


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

void AddTimeDerivative() override;

void CalculateVectorMaterialDerivative() override;

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

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
    StandardRecoveryProcess() = delete;

    /// Assignment operator.
    StandardRecoveryProcess& operator=(StandardRecoveryProcess const& rOther) = delete;

    /// Copy constructor.
    StandardRecoveryProcess(StandardRecoveryProcess const& rOther) = delete;

    ///@}

}; // Class StandardRecoveryProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_STANDARD_RECOVERY_PROCESS_H
