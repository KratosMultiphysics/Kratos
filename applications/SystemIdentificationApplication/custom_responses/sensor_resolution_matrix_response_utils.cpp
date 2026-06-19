//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: SystemIdentificationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/nd_data.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "utilities/data_type_traits.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/filtering/neareset_entity_explicit_damping.h"
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_resolution_matrix_response_utils.h"

namespace Kratos {

SensorResolutionMatrixResponseUtils::SensorResolutionMatrixResponseUtils(
    SensorMaskStatus::Pointer pSensorMaskStatus,
    const double StepSize,
    const double FilterRadius,
    ModelPart& rModelPart,
    const std::string& rKernelFunctionType,
    const IndexType MaxLeafSize,
    const IndexType EchoLevel,
    const bool NodeCloudMesh,
    const bool StoreFilterMatrix)
    : mpSensorMaskStatus(pSensorMaskStatus),
      mStepSize(StepSize)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mStepSize <= 0.0) << "The step size can only be positive value [ step size = " << mStepSize << " ].\n";

    std::visit([this, &rModelPart, &rKernelFunctionType, MaxLeafSize, EchoLevel, NodeCloudMesh, StoreFilterMatrix, FilterRadius](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            auto p_filter = Kratos::make_shared<ExplicitFilterUtils<container_type>>(
                rModelPart, rKernelFunctionType, MaxLeafSize, EchoLevel,
                NodeCloudMesh, StoreFilterMatrix);
            auto p_damper = Kratos::make_shared<NearestEntityExplicitDamping<container_type>>(rModelPart.GetModel(), Parameters("""{}"""), 1);

            auto p_radius_ta = Kratos::make_shared<TensorAdaptor<double>>(pContainer, Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, pContainer->size()), FilterRadius), false);
            p_filter->SetRadius(p_radius_ta);
            p_damper->SetRadius(p_radius_ta);
            p_filter->SetDamping(p_damper);
            p_damper->Update();
            p_filter->Update();
            this->mpFilter = p_filter;
        } else {
            KRATOS_ERROR << "The SensorResolutionMatrixResponseUtils only supports nodal, condition or elemental masks.";
        }
    }, pSensorMaskStatus->pGetMaskContainer());

    KRATOS_CATCH("");
}

double SensorResolutionMatrixResponseUtils::CalculateValue()
{
    KRATOS_TRY

    const auto& r_mask_status = mpSensorMaskStatus->GetMaskStatuses();
    mResolutionMatrix.resize(r_mask_status.size1(), r_mask_status.size1(), false);

    return std::visit([&](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            auto p_filter = std::get<typename ExplicitFilterUtils<container_type>::Pointer>(mpFilter);

            double frobenius_norm = 0.0;
            const double coeff = 1 / mStepSize;

            Vector aux_vec(r_mask_status.size1());
            NDData<double>::Pointer p_nd_data = Kratos::make_shared<NDData<double>>(&aux_vec[0], DenseVector<unsigned int>(1, pContainer->size()), false);
            TensorAdaptor<double> tensor_adaptor(pContainer, p_nd_data, false);

            for (IndexType i_col = 0; i_col < mResolutionMatrix.size1(); ++i_col) {
                IndexPartition<IndexType>(mResolutionMatrix.size1()).for_each([&r_mask_status, &aux_vec, i_col](const auto i_row){
                    double& value = aux_vec[i_row];
                    value = 0.0;
                    for (IndexType k = 0; k < r_mask_status.size2(); ++k) {
                        value += r_mask_status(i_row, k) * r_mask_status(i_col, k);
                    }
                });

                auto p_filtered_tensor_adaptor = p_filter->ForwardFilterField(*p_filter->BackwardFilterField(tensor_adaptor));
                const auto data_view = p_filtered_tensor_adaptor->ViewData();

                frobenius_norm += IndexPartition<IndexType>(data_view.size()).for_each<SumReduction<double>>([&](const auto iRow){
                    const double value = data_view[iRow];
                    mResolutionMatrix(iRow, i_col) = value;
                    return value * value;
                });

                frobenius_norm += coeff * coeff - 2.0 * coeff * data_view[i_col];
            }

            return frobenius_norm * mStepSize * mStepSize * 0.5;
        } else {
            KRATOS_ERROR << "The SensorResolutionMatrixResponseUtils only supports nodal, condition or elemental masks.";
            return 0.0;
        }
    }, mpSensorMaskStatus->pGetMaskContainer());

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorResolutionMatrixResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const auto& r_masks = mpSensorMaskStatus->GetMasks();

    Matrix auxiliary_mask_matrix(r_masks.size1(), r_masks.size1());

    auto p_sensor_nodes = mpSensorMaskStatus->pGetSensorModelPart()->pNodes();
    auto result_nd_data = Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, p_sensor_nodes->size()));
    auto result_data_view = result_nd_data->ViewData();

    std::visit([&](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            auto p_filter = std::get<typename ExplicitFilterUtils<container_type>::Pointer>(mpFilter);

            const double coeff = 1 / mStepSize;

            for (IndexType i_sensor = 0; i_sensor < p_sensor_nodes->size(); ++i_sensor) {
                const Vector& r_sensor_mask = column(r_masks, i_sensor);
                noalias(auxiliary_mask_matrix) = outer_prod(r_sensor_mask, r_sensor_mask);

                auto& value = result_data_view[i_sensor];
                value = 0.0;

                for (IndexType i_col = 0; i_col < mResolutionMatrix.size2(); ++i_col) {
                    double* p_row_start = &auxiliary_mask_matrix(i_col, 0);
                    NDData<double>::Pointer p_nd_data = Kratos::make_shared<NDData<double>>(p_row_start, DenseVector<unsigned int>(1, r_masks.size1()), false);
                    TensorAdaptor<double> tensor_adaptor(pContainer, p_nd_data, false);
                    auto p_filtered_tensor_adaptor = p_filter->ForwardFilterField(*p_filter->BackwardFilterField(tensor_adaptor));
                    const auto data_view = p_filtered_tensor_adaptor->ViewData();

                    value += IndexPartition<IndexType>(mResolutionMatrix.size2()).for_each<SumReduction<double>>([&](const auto iRow) {
                        return mResolutionMatrix(iRow, i_col) * data_view[iRow];
                    });
                    value -= data_view[i_col] * coeff;
                }

                value *= 2.0 * mStepSize * mStepSize * (mpSensorMaskStatus->GetSensorModelPart().NodesBegin() + i_sensor)->GetValue(SENSOR_STATUS);
            }
        } else {
            KRATOS_ERROR << "The SensorResolutionMatrixResponseUtils only supports nodal, condition or elemental masks.";
        }
    }, mpSensorMaskStatus->pGetMaskContainer());

    return Kratos::make_shared<TensorAdaptor<double>>(p_sensor_nodes, result_nd_data, false);

    KRATOS_CATCH("");
}

std::variant<
    ExplicitFilterUtils<ModelPart::NodesContainerType>::Pointer,
    ExplicitFilterUtils<ModelPart::ConditionsContainerType>::Pointer,
    ExplicitFilterUtils<ModelPart::ElementsContainerType>::Pointer> SensorResolutionMatrixResponseUtils::GetFilter()
{
    return mpFilter;
}

} /* namespace Kratos.*/