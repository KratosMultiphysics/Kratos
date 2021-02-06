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
#include <cmath>
#include <tuple>

// External includes

// Project includes
#include "includes/communicator.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "rans_variable_difference_norm_calculation_utility.h"

namespace Kratos
{
template <typename TDataType>
void RansVariableDifferenceNormsCalculationUtility<TDataType>::InitializeCalculation()
{
    KRATOS_TRY

    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_local_mesh = r_communicator.LocalMesh();

    const auto& r_local_elements = r_local_mesh.Elements();
    const int number_of_elements = r_local_elements.size();
    if (static_cast<int>(mElementData.size()) < number_of_elements) {
        mElementData.resize(number_of_elements);
    }

    IndexPartition<int>(number_of_elements).for_each([&](const int iElement) {
        const auto p_element = r_local_elements.begin() + iElement;
        mElementData[iElement] = p_element->GetValue(mrVariable);
    });

    const auto& r_local_conditions = r_local_mesh.Conditions();
    const int number_of_conditions = r_local_conditions.size();
    if (static_cast<int>(mConditionData.size()) < number_of_conditions) {
        mConditionData.resize(number_of_conditions);
    }

    IndexPartition<int>(number_of_conditions).for_each([&](const int iCondition) {
        const auto p_condition = r_local_conditions.begin() + iCondition;
        mConditionData[iCondition] = p_condition->GetValue(mrVariable);
    });

    KRATOS_CATCH("");
}

template <typename TDataType>
std::tuple<double, double> RansVariableDifferenceNormsCalculationUtility<TDataType>::CalculateDifferenceNorm()
{
    KRATOS_TRY

    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_local_mesh = r_communicator.LocalMesh();

    const auto& r_local_elements = r_local_mesh.Elements();
    const int number_of_elements = r_local_elements.size();

    KRATOS_ERROR_IF(static_cast<int>(mElementData.size()) < number_of_elements)
        << "Data is not properly initialized for " << mrVariable.Name() << " in "
        << mrModelPart.Name() << ". Please use \"InitializeCalculation\" first.\n";

    const auto& r_local_conditions = r_local_mesh.Conditions();
    const int number_of_conditions = r_local_conditions.size();

    KRATOS_ERROR_IF(static_cast<int>(mConditionData.size()) < number_of_conditions)
        << "Data is not properly initialized for " << mrVariable.Name() << " in "
        << mrModelPart.Name() << ". Please use \"InitializeCalculation\" first.\n";

    double element_dx, element_solution;
    std::tie(element_dx, element_solution) =
        IndexPartition<int>(number_of_elements)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(
                [&](const int iElement) -> std::tuple<double, double> {
                    const auto& r_element = *(r_local_elements.begin() + iElement);
                    const double value = r_element.GetValue(mrVariable);

                    return std::make_tuple<double, double>(
                        std::pow(value - mElementData[iElement], 2), std::pow(value, 2));
                });

    double condition_dx, condition_solution;
    std::tie(condition_dx, condition_solution) =
        IndexPartition<int>(number_of_conditions)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(
                [&](const int iCondition) -> std::tuple<double, double> {
                    const auto& r_condition = *(r_local_conditions.begin() + iCondition);
                    const double value = r_condition.GetValue(mrVariable);

                    return std::make_tuple<double, double>(
                        std::pow(value - mConditionData[iCondition], 2), std::pow(value, 2));
                });

    const std::vector<double> norm_values = {element_dx + condition_dx, element_solution + condition_solution, static_cast<double>(number_of_elements + number_of_conditions)};
    const auto& total_norm_values = r_communicator.GetDataCommunicator().SumAll(norm_values);

    const double dx = std::sqrt(total_norm_values[0]);
    double solution = std::sqrt(total_norm_values[1]);
    solution = (solution == 0.0 ? 1.0 : solution);

    return std::make_tuple<double, double>(dx / solution, dx / total_norm_values[2]);

    KRATOS_CATCH("");
}

// template instantiations

template class RansVariableDifferenceNormsCalculationUtility<double>;

} // namespace Kratos.
