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

#include "dof_utilities.h"
#include "includes/variables.h"

using namespace Kratos;

namespace
{

const Variable<double>& GetTimeDerivativeOf(const Dof<double>* pDof)
{
    // Member function `GetVariable` of class `Dof<double>` returns a reference to a `VariableData`
    // instance. However, that class does not give us access to the first time derivative.
    // Therefore, we need to downcast to a `Variable<double>` first. In practice, we expect this
    // cast to always succeed, but we will check it here, just in case.
    auto p_variable = dynamic_cast<const Variable<double>*>(&(pDof->GetVariable()));
    KRATOS_ERROR_IF_NOT(p_variable) << "Variable associated with DOF is not of type double" << std::endl;
    return p_variable->GetTimeDerivative();
}

} // namespace

namespace Kratos::Geo::DofUtilities
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs)
{
    std::vector<std::size_t> result;
    std::transform(rDofs.begin(), rDofs.end(), std::back_inserter(result),
                   [](const auto p_dof) { return p_dof->EquationId(); });
    return result;
}

Vector ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result = Vector(rDofs.size());
    std::transform(rDofs.begin(), rDofs.end(), result.begin(),
                   [BufferIndex](auto p_dof) { return p_dof->GetSolutionStepValue(BufferIndex); });
    return result;
}

Vector ExtractFirstTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result                    = Vector(rDofs.size());
    auto get_first_time_derivative = [BufferIndex](const auto p_dof) -> double {
        return p_dof->GetSolutionStepValue(GetTimeDerivativeOf(p_dof), BufferIndex);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_first_time_derivative);
    return result;
}

Vector ExtractSecondTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result                     = Vector(rDofs.size());
    auto get_second_time_derivative = [BufferIndex](const auto p_dof) -> double {
        return p_dof->GetSolutionStepValue(GetTimeDerivativeOf(p_dof).GetTimeDerivative(), BufferIndex);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_second_time_derivative);
    return result;
}

Vector ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result        = Vector(rDofs.size());
    auto get_dof_value = [BufferIndex](const auto p_dof) -> double {
        // Why should we return 0.0 for water pressures?
        return p_dof->GetVariable() == WATER_PRESSURE ? 0.0 : p_dof->GetSolutionStepValue(BufferIndex);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_dof_value);
    return result;
}

Vector ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result                    = Vector(rDofs.size());
    auto get_first_time_derivative = [BufferIndex](const auto p_dof) -> double {
        // Why should we return 0.0 for the first time derivatives of water pressure?
        if (p_dof->GetVariable() == WATER_PRESSURE) return 0.0;

        return p_dof->GetSolutionStepValue(GetTimeDerivativeOf(p_dof), BufferIndex);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_first_time_derivative);
    return result;
}

Vector ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex)
{
    auto result                     = Vector(rDofs.size());
    auto get_second_time_derivative = [BufferIndex](const auto p_dof) -> double {
        // Why should we return 0.0 for the second time derivatives of water pressure?
        if (p_dof->GetVariable() == WATER_PRESSURE) return 0.0;

        return p_dof->GetSolutionStepValue(GetTimeDerivativeOf(p_dof).GetTimeDerivative(), BufferIndex);
    };
    std::transform(rDofs.begin(), rDofs.end(), result.begin(), get_second_time_derivative);
    return result;
}

} // namespace Kratos::Geo::DofUtilities
