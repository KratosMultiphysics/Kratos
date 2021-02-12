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
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "rans_variable_difference_norm_calculation_utility.h"

namespace Kratos
{
template <typename TDataType>
void RansVariableDifferenceNormsCalculationUtility<TDataType>::InitializeCalculation()
{
    KRATOS_TRY

    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const auto& r_local_mesh = r_communicator.LocalMesh();

    const auto& r_local_elements = r_local_mesh.Elements();
    const int number_of_elements = r_local_elements.size();
    if (static_cast<int>(mElementData.size()) < number_of_elements) {
        mElementData.resize(number_of_elements);
    }

    using tls_type = std::tuple<Vector, Matrix, Geometry<Node<3>>::ShapeFunctionsGradientsType>;
    IndexPartition<int>(number_of_elements).for_each(tls_type(), [&](const int iElement, tls_type& rTLS) {
        const auto& r_element = *(r_local_elements.begin() + iElement);

        auto& Ws = std::get<0>(rTLS);
        auto& Ns = std::get<1>(rTLS);
        auto& dNdXs = std::get<2>(rTLS);

        RansCalculationUtilities::CalculateGeometryData(r_element.GetGeometry(), GeometryData::IntegrationMethod::GI_GAUSS_1, Ws, Ns, dNdXs);

        ConstitutiveLaw::Parameters parameters(
            r_element.GetGeometry(), r_element.GetProperties(), r_process_info);
        parameters.SetShapeFunctionsValues(row(Ns, 0));
        parameters.SetShapeFunctionsDerivatives(dNdXs[0]);

        auto p_constitutive_law = r_element.GetValue(CONSTITUTIVE_LAW);

        p_constitutive_law->CalculateValue(parameters, mrVariable, mElementData[iElement]);
    });

    KRATOS_CATCH("");
}

template <typename TDataType>
std::tuple<double, double> RansVariableDifferenceNormsCalculationUtility<TDataType>::CalculateDifferenceNorm()
{
    KRATOS_TRY

    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const auto& r_local_mesh = r_communicator.LocalMesh();

    const auto& r_local_elements = r_local_mesh.Elements();
    const int number_of_elements = r_local_elements.size();

    KRATOS_ERROR_IF(static_cast<int>(mElementData.size()) < number_of_elements)
        << "Data is not properly initialized for " << mrVariable.Name() << " in "
        << mrModelPart.Name() << ". Please use \"InitializeCalculation\" first.\n";

    using tls_type = std::tuple<Vector, Matrix, Geometry<Node<3>>::ShapeFunctionsGradientsType>;

    double element_dx, element_solution;
    std::tie(element_dx, element_solution) =
        IndexPartition<int>(number_of_elements)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(tls_type(),
                [&](const int iElement, tls_type& rTLS) -> std::tuple<double, double> {
                    const auto& r_element = *(r_local_elements.begin() + iElement);

                    auto& Ws = std::get<0>(rTLS);
                    auto& Ns = std::get<1>(rTLS);
                    auto& dNdXs = std::get<2>(rTLS);

                    RansCalculationUtilities::CalculateGeometryData(r_element.GetGeometry(), GeometryData::IntegrationMethod::GI_GAUSS_1, Ws, Ns, dNdXs);

                    ConstitutiveLaw::Parameters parameters(
                        r_element.GetGeometry(), r_element.GetProperties(), r_process_info);
                    parameters.SetShapeFunctionsValues(row(Ns, 0));
                    parameters.SetShapeFunctionsDerivatives(dNdXs[0]);

                    auto p_constitutive_law = r_element.GetValue(CONSTITUTIVE_LAW);

                    double value;
                    p_constitutive_law->CalculateValue(parameters, mrVariable, value);

                    return std::make_tuple<double, double>(
                        std::pow(value - mElementData[iElement], 2), std::pow(value, 2));
                });

    const std::vector<double> norm_values = {element_dx, element_solution, static_cast<double>(number_of_elements)};
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
