// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include <algorithm>

#include "dof_utilities.h"
#include "includes/variables.h"

namespace Kratos
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs)
{
    std::vector<std::size_t> result;
    std::transform(rDofs.begin(), rDofs.end(), std::back_inserter(result),
                   [](const auto p_dof) { return p_dof->EquationId(); });
    return result;
}

Vector ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int Step)
{
    auto result = Vector(rDofs.size());
    std::transform(rDofs.begin(), rDofs.end(), result.begin(),
                   [Step](auto p_dof) { return p_dof->GetSolutionStepValue(Step); });
    return result;
}

Vector ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step)
{
    auto result        = Vector(rDofs.size());
    auto get_dof_value = [Step](const auto p_dof) -> double {
        // Why should we return 0.0 for water pressures?
        return p_dof->GetVariable() == WATER_PRESSURE ? 0.0 : p_dof->GetSolutionStepValue(Step);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_dof_value);
    return result;
}

Vector ExtractFirstDerivatives(const std::vector<Dof<double>*>& rDofs, int Step)
{
    auto result                    = Vector(rDofs.size());
    auto get_first_time_derivative = [Step](const auto p_dof) -> double {
        // Unfortunately, the first time derivative cannot be accessed from a `VariableData`
        // instance. However, we assume that each DOF corresponds to a `Variable<double>` instance.
        // Then we should be able to downcast the `VariableData` instance.
        auto p_variable = dynamic_cast<const Variable<double>*>(&p_dof->GetVariable());
        KRATOS_ERROR_IF_NOT(p_variable) << "Variable associated with DOF is not of type double" << std::endl;
        return p_dof->GetSolutionStepValue(p_variable->GetTimeDerivative(), Step);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_first_time_derivative);
    return result;
}

Vector ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step)
{
    auto result                    = Vector(rDofs.size());
    auto get_first_time_derivative = [Step](const auto p_dof) -> double {
        // Why should we return 0.0 for the first time derivatives of water pressure?
        if (p_dof->GetVariable() == WATER_PRESSURE) return 0.0;

        // Unfortunately, the first time derivative cannot be accessed from a `VariableData`
        // instance. However, we know that for all U-Pw elements the degrees of freedom are either
        // displacement components or water pressures, which correspond to `Variable<double>`
        // instances. Therefore, we should be able to downcast the `VariableData` instance.
        auto p_variable = dynamic_cast<const Variable<double>*>(&p_dof->GetVariable());
        KRATOS_ERROR_IF_NOT(p_variable) << "Variable associated with DOF is not of type double" << std::endl;
        return p_dof->GetSolutionStepValue(p_variable->GetTimeDerivative(), Step);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_first_time_derivative);
    return result;
}

Vector ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int Step)
{
    auto result                     = Vector(rDofs.size());
    auto get_second_time_derivative = [Step](const auto p_dof) -> double {
        // Why should we return 0.0 for the second time derivatives of water pressure?
        if (p_dof->GetVariable() == WATER_PRESSURE) return 0.0;

        // Unfortunately, the first time derivative cannot be accessed from a `VariableData`
        // instance. However, we know that for all U-Pw elements the degrees of freedom are either
        // displacement components or water pressures, which correspond to `Variable<double>`
        // instances. Therefore, we should be able to downcast the `VariableData` instance.
        auto p_variable = dynamic_cast<const Variable<double>*>(&p_dof->GetVariable());
        KRATOS_ERROR_IF_NOT(p_variable) << "Variable associated with DOF is not of type double" << std::endl;
        return p_dof->GetSolutionStepValue(p_variable->GetTimeDerivative().GetTimeDerivative(), Step);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_second_time_derivative);
    return result;
}

} // namespace Kratos
