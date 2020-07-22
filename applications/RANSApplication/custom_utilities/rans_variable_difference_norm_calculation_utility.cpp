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
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();
    const int number_of_nodes = r_nodes.size();

    KRATOS_ERROR_IF(!mrModelPart.HasNodalSolutionStepVariable(mrVariable))
        << mrVariable.Name() << " not found in nodal solution step variables list in "
        << mrModelPart.Name() << ".\n";

    if (static_cast<int>(mData.size()) < number_of_nodes) {
        mData.resize(number_of_nodes);
    }

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = *(r_nodes.begin() + i_node);
        mData[i_node] = r_node.FastGetSolutionStepValue(mrVariable);
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
std::tuple<double, double> RansVariableDifferenceNormsCalculationUtility<TDataType>::CalculateDifferenceNorm()
{
    KRATOS_TRY

    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();
    const int number_of_nodes = r_nodes.size();

    KRATOS_ERROR_IF(static_cast<int>(mData.size()) < number_of_nodes)
        << "Data is not properly initialized for " << mrVariable.Name() << " in "
        << mrModelPart.Name() << ". Please use \"InitializeCalculation\" first.\n";

    double dx{0.0}, solution{0.0};
#pragma omp parallel for reduction(+ : dx, solution)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = *(r_nodes.begin() + i_node);
        const double value = r_node.FastGetSolutionStepValue(mrVariable);
        dx += std::pow(value - mData[i_node], 2);
        solution += std::pow(value, 2);
    }

    const std::vector<double> norm_values = {
        dx, solution, static_cast<double>(number_of_nodes)};
    const std::vector<double>& total_norm_values =
        r_communicator.GetDataCommunicator().SumAll(norm_values);

    dx = std::sqrt(total_norm_values[0]);
    solution = std::sqrt(total_norm_values[1]);
    solution = (solution == 0.0 ? 1.0 : solution);

    return std::make_tuple<double, double>(dx / solution, dx / total_norm_values[2]);

    KRATOS_CATCH("");
}

// template instantiations

template class RansVariableDifferenceNormsCalculationUtility<double>;

} // namespace Kratos.
