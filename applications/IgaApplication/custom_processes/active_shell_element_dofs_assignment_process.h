//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class ActiveShellElementDofAssignmentProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) ActiveShellElementDofAssignmentProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ActiveShellElementDofAssignmentProcess
    KRATOS_CLASS_POINTER_DEFINITION(ActiveShellElementDofAssignmentProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ActiveShellElementDofAssignmentProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~ActiveShellElementDofAssignmentProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    /// Called once before the solution loop and is writing the quadrature domain.
    void ExecuteInitialize() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ActiveShellElementDofAssignmentProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ActiveShellElementDofAssignmentProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    /// Model part and different settings

    Model& mrModel;             /// The main model part
    
    std::string mIgaModelPartName;

    std::string mActiveShellDofModelPartName;

    ///@}

}; // Class ActiveShellElementDofAssignmentProcess

///@}
}  // namespace Kratos.

