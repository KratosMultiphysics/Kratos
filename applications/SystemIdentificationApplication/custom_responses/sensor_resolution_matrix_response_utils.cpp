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

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_resolution_matrix_response_utils.h"

namespace Kratos {

SensorResolutionMatrixResponseUtils::SensorResolutionMatrixResponseUtils(
    SensorMaskStatus::Pointer pSensorMaskStatus,
    const double StepSize,
    const double FilterRadius,
    const ModelPart& rModelPart,
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
            p_filter->SetRadius(Kratos::make_shared<TensorAdaptor<double>>(pContainer, Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, pContainer->size()), FilterRadius), false));
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

    Matrix auxiliary_mask_matrix(r_mask_status.size1(), r_mask_status.size1());
    mResolutionMatrix.resize(r_mask_status.size1(), r_mask_status.size1(), false);

    noalias(auxiliary_mask_matrix) = prod(r_mask_status, trans(r_mask_status));

    return std::visit([&](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(IsInList<container_type, ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            auto p_filter = std::get<typename ExplicitFilterUtils<container_type>::Pointer>(mpFilter);

            double frobenius_norm = 0.0;
            const double coeff = 1 / mStepSize;
            for (IndexType i_col = 0; i_col < mResolutionMatrix.size1(); ++i_col) {
                double* p_row_start = &auxiliary_mask_matrix(i_col, 0);
                NDData<double>::Pointer p_nd_data = Kratos::make_shared<NDData<double>>(p_row_start, DenseVector<unsigned int>(1, pContainer->size()), false);
                TensorAdaptor<double> tensor_adaptor(pContainer, p_nd_data, false);
                auto p_filtered_tensor_adaptor = p_filter->ForwardFilterField(*p_filter->BackwardFilterField(tensor_adaptor));
                const auto data_view = p_filtered_tensor_adaptor->ViewData();

                IndexPartition<IndexType>(data_view.size()).for_each([&](const auto iRow){
                    const double value = data_view[iRow];
                    mResolutionMatrix(iRow, i_col) = value;
                    frobenius_norm += value * value;
                });
                const double value = data_view[i_col];
                frobenius_norm = frobenius_norm - (value * value) + (value - coeff) * (value - coeff);
            }

            return frobenius_norm * mStepSize * mStepSize;
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

                for (IndexType i_col = 0; i_col < mResolutionMatrix.size1(); ++i_col) {
                    double* p_row_start = &auxiliary_mask_matrix(i_col, 0);
                    NDData<double>::Pointer p_nd_data = Kratos::make_shared<NDData<double>>(p_row_start, DenseVector<unsigned int>(1, r_masks.size1()), false);
                    TensorAdaptor<double> tensor_adaptor(pContainer, p_nd_data, false);
                    auto p_filtered_tensor_adaptor = p_filter->ForwardFilterField(*p_filter->BackwardFilterField(tensor_adaptor));
                    const auto data_view = p_filtered_tensor_adaptor->ViewData();

                    IndexPartition<IndexType>(mResolutionMatrix.size2()).for_each([&](const auto Index) {
                        value += mResolutionMatrix(i_col, Index) * data_view[Index];
                    });

                    value -= data_view[i_col] * coeff;
                    value *= 2.0 * mStepSize * mStepSize * (mpSensorMaskStatus->GetSensorModelPart().NodesBegin() + i_sensor)->GetValue(SENSOR_STATUS);
                }
            }
        } else {
            KRATOS_ERROR << "The SensorResolutionMatrixResponseUtils only supports nodal, condition or elemental masks.";
        }
    }, mpSensorMaskStatus->pGetMaskContainer());

    return Kratos::make_shared<TensorAdaptor<double>>(p_sensor_nodes, result_nd_data, false);

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/