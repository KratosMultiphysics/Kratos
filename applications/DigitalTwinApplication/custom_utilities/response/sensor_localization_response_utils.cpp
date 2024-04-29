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
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_localization_response_utils.h"

namespace Kratos {

SensorLocalizationResponseUtils::SensorLocalizationResponseUtils(
    ModelPart& rSensorModelPart,
    const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList,
    const double P)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P),
      mMasksList(rMasksList)
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    KRATOS_ERROR_IF(mMasksList.empty())
        << "Masks list is empty.";

    KRATOS_ERROR_IF_NOT(r_model_part.NumberOfNodes() == mMasksList.size())
        << "Number of sensors and masks list mismatch [ number of sensors = "
        << r_model_part.NumberOfNodes() << ", number of masks  = "
        << mMasksList.size() << " ].\n";

    // get the total domain size and domain size ratios
    const auto& r_container = mMasksList.front()->GetContainer();
    mDomainSizeRatio.resize(r_container.size());
    const double local_domain_size = IndexPartition<IndexType>(r_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
        const double domain_size = (r_container.begin() + Index)->GetGeometry().DomainSize();
        mDomainSizeRatio[Index] = domain_size;
        return domain_size;
    });
    const double total_domain_size = mMasksList.front()->GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(local_domain_size);
    block_for_each(mDomainSizeRatio, [total_domain_size](auto& rValue) {
        rValue /= total_domain_size;
    });

    KRATOS_CATCH("");
}

void SensorLocalizationResponseUtils::ComputeClusterDifference(Matrix& rOutput) const
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    const IndexType number_of_elements = mMasksList.front()->GetContainer().size();

    if (rOutput.size1() != number_of_elements || rOutput.size2() != number_of_elements) {
        rOutput.resize(number_of_elements, number_of_elements, false);
    }

    noalias(rOutput) = ZeroMatrix(number_of_elements, number_of_elements);

    IndexPartition<IndexType>(number_of_elements).for_each([&](const auto i) {
        for (IndexType k_sensor = 0; k_sensor < r_model_part.NumberOfNodes(); ++k_sensor) {
            const double sensor_status = (r_model_part.NodesBegin() + k_sensor)->GetValue(SENSOR_STATUS);
            const auto& r_mask_exp = mMasksList[k_sensor]->GetExpression();
            const bool i_value = static_cast<bool>(r_mask_exp.Evaluate(i, i, 0));
            for (IndexType j = i + 1; j < number_of_elements; ++j) {
                const double value = sensor_status * (i_value ^ static_cast<bool>(r_mask_exp.Evaluate(j, j ,0)));
                rOutput(i, j) += value;
                rOutput(j, i) += value;
            }
        }
    });

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue() const
{
    KRATOS_TRY

    SmoothClamper<ModelPart::ElementsContainerType> clamper(0.0, 1.0);

    Matrix aux_matrix;
    ComputeClusterDifference(aux_matrix);

    IndexPartition<IndexType>(aux_matrix.size1()).for_each([&](const auto i) {
        for (IndexType j = i + 1; j <  aux_matrix.size2(); ++j) {
            aux_matrix(i, j) = (1.0 - clamper.Clamp(aux_matrix(i, j)));
            aux_matrix(j, i) = aux_matrix(i, j);
        }
        aux_matrix(i, i) = 1.0;
    });

    const auto& r_mask_container = mMasksList.front()->GetContainer();

    const double summation = IndexPartition<IndexType>(r_mask_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
        double cluster_size_ratio = 0.0;
        const Vector& r_values = row(aux_matrix, Index);
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            cluster_size_ratio += mDomainSizeRatio[i_element] * r_values[i_element];
        }
        return std::pow(cluster_size_ratio, mP);
    });

    return std::pow(summation, 1 / mP);

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorLocalizationResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    auto& r_model_part = *mpSensorModelPart;

    const auto& r_mask_container = mMasksList.front()->GetContainer();

    SmoothClamper<ModelPart::ElementsContainerType> clamper(0.0, 1.0);

    // first compute the cluster differences
    Matrix aux_matrix;
    ComputeClusterDifference(aux_matrix);

    Matrix aux_matrix_2(aux_matrix.size1(), aux_matrix.size2());
    IndexPartition<IndexType>(aux_matrix.size1()).for_each([&](const auto i) {
        for (IndexType j = i + 1; j <  aux_matrix.size2(); ++j) {
            aux_matrix_2(i, j) = (1.0 - clamper.Clamp(aux_matrix(i, j)));
            aux_matrix_2(j, i) = aux_matrix_2(i, j);
        }
        aux_matrix_2(i, i) = 1.0;
    });

    std::vector<double> cluster_size_ratios(r_mask_container.size(), 0.0);
    const double summation = IndexPartition<IndexType>(r_mask_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
        const Vector& r_values = row(aux_matrix_2, Index);
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            cluster_size_ratios[Index] += mDomainSizeRatio[i_element] * r_values[i_element];
        }
        return std::pow(cluster_size_ratios[Index], mP);
    });

    // now we compute the clamper derivatives
    IndexPartition<IndexType>(aux_matrix.size1()).for_each([&](const auto i) {
        for (IndexType j = i + 1; j <  aux_matrix.size2(); ++j) {
            aux_matrix(i, j) = clamper.ClampDerivative(aux_matrix(i, j));
        }
    });

    const double coeff = std::pow(summation, 1 / mP - 1);

    auto p_expression = LiteralFlatExpression<double>::Create(r_model_part.NumberOfNodes(), {});
    auto& r_expression = *p_expression;
    IndexPartition<IndexType>(r_model_part.NumberOfNodes()).for_each([&](const auto SensorIndex) {
        const auto& r_mask_exp = mMasksList[SensorIndex]->GetExpression();
        double& value = *(r_expression.begin() + SensorIndex);
        value = 0.0;
        for (IndexType i_element = 0; i_element < r_mask_container.size(); ++i_element) {
            double d_cluster_size_ratio_d_sensor_status = 0.0;
            const double domain_size_ratio = mDomainSizeRatio[i_element];
            const bool i_value = static_cast<bool>(r_mask_exp.Evaluate(i_element, i_element, 0));

            const Vector& aux_vec_1 = column(aux_matrix, i_element);
            for (IndexType j_element = 0; j_element < i_element; ++j_element) {
                d_cluster_size_ratio_d_sensor_status -= domain_size_ratio * aux_vec_1[j_element] * (i_value ^ static_cast<bool>(r_mask_exp.Evaluate(j_element, j_element, 0)));
            }

            const Vector& aux_vec_2 = row(aux_matrix, i_element);
            for (IndexType j_element = i_element + 1; j_element < r_mask_container.size(); ++j_element) {
                d_cluster_size_ratio_d_sensor_status -= domain_size_ratio * aux_vec_2[j_element] * (i_value ^ static_cast<bool>(r_mask_exp.Evaluate(j_element, j_element, 0)));
            }

            value += std::pow(cluster_size_ratios[i_element], mP - 1) * d_cluster_size_ratio_d_sensor_status;
        }
        value *= coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(r_model_part);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/