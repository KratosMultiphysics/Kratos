//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_localization_response_utils.h"

namespace Kratos {

SensorLocalizationResponseUtils::SensorLocalizationResponseUtils(
    ModelPart& rSensorModelPart,
    const double P)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P)
{
}


void SensorLocalizationResponseUtils::ComputeClusterDifference(
    Matrix& rOutput,
    const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList)
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    const IndexType number_of_elements = rMasksList.front()->GetContainer().size();

    if (rOutput.size1() != number_of_elements || rOutput.size2() != number_of_elements) {
        rOutput.resize(number_of_elements, number_of_elements, false);
    }

    noalias(rOutput) = ZeroMatrix(number_of_elements, number_of_elements);

    IndexPartition<IndexType>(number_of_elements).for_each([&rOutput, &rMasksList, &r_model_part, number_of_elements](const auto i) {
        for (IndexType i_sensor = 0; i_sensor < r_model_part.NumberOfNodes(); ++i_sensor) {
            const auto& r_mask_exp = rMasksList[i_sensor]->GetExpression();
            const double sensor_status = (r_model_part.NodesBegin() + i_sensor)->GetValue(SENSOR_STATUS);
            const double i_value = r_mask_exp.Evaluate(i, i, 0);
            for (IndexType j = i + 1; j < number_of_elements; ++j) {
                const double value = sensor_status * std::pow(i_value - r_mask_exp.Evaluate(j, j, 0), 2);
                rOutput(i, j) += value;
                rOutput(j, i) += value;
            }
        }
    });

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue(const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList)
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    KRATOS_ERROR_IF(rMasksList.empty())
        << "Masks list is empty.";

    KRATOS_ERROR_IF_NOT(r_model_part.NumberOfNodes() == rMasksList.size())
        << "Number of sensors and masks list mismatch [ number of sensors = "
        << r_model_part.NumberOfNodes() << ", number of masks  = "
        << rMasksList.size() << " ].\n";

    SmoothClamper<ModelPart::ElementsContainerType> clamper(0.0, 1.0);

    Matrix cluster_differences;
    ComputeClusterDifference(cluster_differences, rMasksList);

    const auto& r_mask_container = rMasksList.front()->GetContainer();

    double summation, total_domain_size;
    std::tie(summation, total_domain_size) = IndexPartition<IndexType>(r_mask_container.size()).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto Index) {
        double cluster_size = 0.0;
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            const double domain_size = (r_mask_container.begin() + i_element)->GetGeometry().DomainSize();
            cluster_size += domain_size * (1.0 - clamper.Clamp(cluster_differences(Index, i_element)));
        }
        return std::make_tuple(std::pow(cluster_size, mP), (r_mask_container.begin() + Index)->GetGeometry().DomainSize());
    });

    return std::pow(summation, 1 / mP) / total_domain_size;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorLocalizationResponseUtils::CalculateGradient(const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList)
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    KRATOS_ERROR_IF(rMasksList.empty())
        << "Masks list is empty.";

    const auto& r_mask_container = rMasksList.front()->GetContainer();

    KRATOS_ERROR_IF_NOT(r_model_part.NumberOfNodes() == rMasksList.size())
        << "Number of sensors and masks list mismatch [ number of sensors = "
        << r_model_part.NumberOfNodes() << ", number of masks  = "
        << rMasksList.size() << " ].\n";

    SmoothClamper<ModelPart::ElementsContainerType> clamper(0.0, 1.0);

    // first compute the cluster differences
    Matrix cluster_differences;
    ComputeClusterDifference(cluster_differences, rMasksList);

    // then compute the cluster sizes
    std::vector<double> cluster_sizes(r_mask_container.size(), 0.0);
    double summation, total_domain_size;
    std::tie(summation, total_domain_size) = IndexPartition<IndexType>(r_mask_container.size()).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto Index) {
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            const double domain_size = (r_mask_container.begin() + i_element)->GetGeometry().DomainSize();
            cluster_sizes[Index] += domain_size * (1.0 - clamper.Clamp(cluster_differences(Index, i_element)));
        }
        return std::make_tuple(std::pow(cluster_sizes[Index], mP), (r_mask_container.begin() + Index)->GetGeometry().DomainSize());
    });

    const double coeff = std::pow(summation, 1 / mP - 1) / (mP * total_domain_size);;

    auto p_expression = LiteralFlatExpression<double>::Create(r_model_part.NumberOfNodes(), {});
    auto& r_expression = *p_expression;
    IndexPartition<IndexType>(r_model_part.NumberOfNodes()).for_each([&](const auto Index) {
        const auto& r_mask_exp = rMasksList[Index]->GetExpression();
        double& value = *(r_expression.begin() + Index);
        value = 0.0;
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            const double i_value = r_mask_exp.Evaluate(i_element, i_element, 0);
            double d_cluster_size_d_sensor_status = 0.0;
            for (IndexType j_element = 0; j_element < r_mask_container.size(); ++j_element) {
                const double domain_size = (r_mask_container.begin() + i_element)->GetGeometry().DomainSize();
                d_cluster_size_d_sensor_status -= domain_size * clamper.ClampDerivative(cluster_differences(i_element, j_element)) * std::pow(i_value - r_mask_exp.Evaluate(j_element, j_element, 0), 2);
            }
            value += mP * std::pow(cluster_sizes[i_element], mP - 1) * d_cluster_size_d_sensor_status;
        }

        value *= coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(r_model_part);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/