//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//

#if !defined(KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Computes the reaction forces for slip modelpart, can be
 *        further used for drag calculation
 *
 * This process sets epsilon values based on the following formula
 *
 * \[
 *
 *  \REACTION = p.n + tau
 *
 * \]
 *
 *
 */

class KRATOS_API(RANS_APPLICATION) RansComputeReactionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansComputeReactionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansComputeReactionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansComputeReactionsProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansComputeReactionsProcess() override = default;

    /// Assignment operator.
    RansComputeReactionsProcess& operator=(
        RansComputeReactionsProcess const& rOther) = delete;

    /// Copy constructor.
    RansComputeReactionsProcess(RansComputeReactionsProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    int Check() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;

    int mEchoLevel;
    bool mPeriodic;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculates reaction on condition nodes
     *
     * Calculates reaction based on condition FRICTION_VELOCITY and PRESSURE values
     * where wall laws are activated.
     *
     * @param rCondition
     */
    void CalculateReactionValues(
        ModelPart::ConditionType& rCondition);

    /**
     * @brief Corrects given variable for periodic conditions
     *
     * This method adds contributions from its periodic neighbour node's rVariable value
     * to own rVariable value and averages them.
     *
     * @param rModelPart
     * @param rVariable
     */
    void CorrectPeriodicNodes(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable);

    ///@}

}; // Class RansComputeReactionsProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansComputeReactionsProcess& rThis);

///@}
///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_COMPUTE_REACTIONS_PROCESS_H_INCLUDED defined