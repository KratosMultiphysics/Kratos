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

// System includes
#include <sstream>
#include <string>
#include <vector>

// External includes

// Project includes

// Include base h
#include "rans_point_execution_formulation_process.h"

namespace Kratos
{
void RansPointExecutionFormulationProcess::ExecuteInitialize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::INITIALIZE_SOLUTION_STEP) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::ExecuteBeforeCouplingSolveStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::BEFORE_COUPLING_SOLVE_STEP) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::Execute()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::EXECUTE) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::ExecuteAfterCouplingSolveStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::AFTER_COUPLING_SOLVE_STEP) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::ExecuteFinalizeSolutionStep()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE_SOLUTION_STEP) {
            this->ExecuteOperation();
            break;
        }
    }
}

void RansPointExecutionFormulationProcess::ExecuteFinalize()
{
    for (const auto& execution_point : mExecutionPointsList) {
        if (execution_point == ExecutionPoint::FINALIZE) {
            this->ExecuteOperation();
            break;
        }
    }
}

std::string RansPointExecutionFormulationProcess::Info() const
{
    return std::string("RansPointExecutionFormulationProcess");
}

void RansPointExecutionFormulationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansPointExecutionFormulationProcess::PrintData(std::ostream& rOStream) const
{
}

void RansPointExecutionFormulationProcess::UpdateExecutionPointsList(
    const std::vector<std::string>& rCopyExecutionPointsList)
{
    KRATOS_TRY

    mExecutionPointsList.clear();

    for (const auto& copy_execution_point : rCopyExecutionPointsList) {
        if (copy_execution_point == "initialize") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE);
        } else if (copy_execution_point == "initialize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::INITIALIZE_SOLUTION_STEP);
        } else if (copy_execution_point == "before_coupling_solve_step") {
            mExecutionPointsList.push_back(ExecutionPoint::BEFORE_COUPLING_SOLVE_STEP);
        } else if (copy_execution_point == "execute") {
            mExecutionPointsList.push_back(ExecutionPoint::EXECUTE);
        } else if (copy_execution_point == "after_coupling_solve_step") {
            mExecutionPointsList.push_back(ExecutionPoint::AFTER_COUPLING_SOLVE_STEP);
        } else if (copy_execution_point == "finalize_solution_step") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE_SOLUTION_STEP);
        } else if (copy_execution_point == "finalize") {
            mExecutionPointsList.push_back(ExecutionPoint::FINALIZE);
        } else {
            KRATOS_ERROR << "Unsupported copy execution point provided. [ "
                            "copy_execution_point = "
                         << copy_execution_point << " ]. Supported points are: \n\n"
                         << "\tinitialize\n"
                         << "\tinitialize_solution_step\n"
                         << "\tbefore_coupling_solve_step\n"
                         << "\texecute\n"
                         << "\tafter_coupling_solve_step\n"
                         << "\tfinalize_solution_step\n"
                         << "\tfinalize\n";
        }
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
