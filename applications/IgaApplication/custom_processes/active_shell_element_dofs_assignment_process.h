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

/**
 * @class ActiveShellElementDofAssignmentProcess
 * @ingroup IgaApplication
 * @brief Creates and configures actuation DOFs for the active shell element from input.
 * @details This process creates a dedicated model part containing one "active shell node" per element
 *          parent geometry (stored as a GlobalPointer in the geometry under the variable
 *          ACTIVE_SHELL_NODE_GP). The node holds the actuation variables (alpha, beta, gamma,
 *          kappa_1, kappa_2, kappa_12) and their corresponding reaction/adjoint variables.
 *
 *          The actuation values and the fix/free status are controlled by three parallel lists:
 *          - `applied_actuation_list`: list of keys ("alpha", "beta", "gamma", "kappa_1", "kappa_2", "kappa_12")
 *          - `applied_actuation_value`: list of values (same length as `applied_actuation_list`)
 *          - `unfixed_actuation_list`: list of flags ("fix" or "free"), controlling whether the DOF is fixed
 *            (Dirichlet) or left free.
 *
 *          @par Example (ProjectParameters.json)
 *          @code
 *              {
 *               "iga_model_part_name": "IgaModelPart",
 *               "active_shell_model_part_name": "ActiveShellDofs",
 *               "applied_actuation_list": ["alpha", "kappa_12"],
 *               "applied_actuation_value": [0.5, 0.1],
 *               "unfixed_actuation_list": ["fix", "free"]
 *              }
 *          @endcode
 */
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
    /**
     * @brief Construct the process from input parameters.
     * @param rModel The global model holding all model parts.
     * @param ThisParameters Settings of this process.
     * @details Expected keys:
     *          - `iga_model_part_name`: model part containing the shell elements to be actuated.
     *          - `active_shell_model_part_name`: name of the new model part that will store the created actuation nodes.
     *          - `applied_actuation_list`: list of actuation keys.
     *          - `applied_actuation_value`: list of values (same length as `applied_actuation_list`).
     *          - `unfixed_actuation_list`: list of "fix"/"free" flags (same length as `applied_actuation_list`).
     */
    ActiveShellElementDofAssignmentProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~ActiveShellElementDofAssignmentProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    /// Called before initialize solution loop, creates the active shell model part with one global node, containing variables, and actuation DOFs.
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

    std::vector<std::string> mAppliedActuationList;

    Vector mAppliedActuationValue;

    std::vector<std::string> mUnfixedActuationList;

    ///@}

}; // Class ActiveShellElementDofAssignmentProcess

///@}
}  // namespace Kratos.

