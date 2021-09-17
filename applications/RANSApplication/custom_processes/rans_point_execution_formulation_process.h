//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_POINT_EXECUTION_FORMULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_POINT_EXECUTION_FORMULATION_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_processes/rans_formulation_process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This class is extending standard Process interface
 *
 * This class extends standard Process interface to have methods to be
 * called before and after coupling step solve in Formulations. This is inherited
 * from Process so, same processes can be used as normal Processes, or if required
 * they can be used in Formulation Process lists.
 *
 */
class RansPointExecutionFormulationProcess : public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansPointExecutionFormulationProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansPointExecutionFormulationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RansPointExecutionFormulationProcess() : RansFormulationProcess()
    {
    }

    /// Destructor.
    ~RansPointExecutionFormulationProcess() override = default;

    RansPointExecutionFormulationProcess& operator=(
        RansPointExecutionFormulationProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteBeforeCouplingSolveStep() override;

    void Execute() override;

    void ExecuteAfterCouplingSolveStep() override;

    void ExecuteFinalizeSolutionStep() override;

    void ExecuteFinalize() override;

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

protected:
    ///@name Protected Member Variables
    ///@{

    enum ExecutionPoint
    {
        INITIALIZE = 0,
        INITIALIZE_SOLUTION_STEP = 1,
        BEFORE_COUPLING_SOLVE_STEP = 2,
        EXECUTE = 3,
        AFTER_COUPLING_SOLVE_STEP = 4,
        FINALIZE_SOLUTION_STEP = 5,
        FINALIZE = 6
    };

    std::vector<ExecutionPoint> mExecutionPointsList;

    ///@}
    ///@name Protected Operations
    ///@{

    void UpdateExecutionPointsList(const std::vector<std::string>& rExecutionPointsList);

    virtual void ExecuteOperation() {};

    ///@}

}; // Class RansPointExecutionFormulationProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(
    std::istream& rIStream,
    RansPointExecutionFormulationProcess& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansPointExecutionFormulationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_RANS_POINT_EXECUTION_FORMULATION_PROCESS_H_INCLUDED  defined
