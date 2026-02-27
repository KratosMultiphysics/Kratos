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

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_inverse_distance_summation_response_utils.h"

namespace Kratos
{

double SensorInverseDistanceSummationResponseUtils::CalculateValue(
    ModelPart &rModelPart,
    const double P,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    return IndexPartition<IndexType>(rDistanceMatrix.GetEntriesSize()).for_each<SumReduction<double>>([P, &rModelPart, &rDistanceMatrix](const auto Index) {
        const auto& index_pair = rDistanceMatrix.GetIndexPair(Index);

        const auto i_index = std::get<0>(index_pair);
        const auto j_index = std::get<1>(index_pair);

        if (i_index != j_index) {
            return std::pow((rModelPart.NodesBegin() + i_index)->GetValue(SENSOR_STATUS), P) * std::pow((rModelPart.NodesBegin() + j_index)->GetValue(SENSOR_STATUS), P) / rDistanceMatrix.GetDistance(Index);
        } else {
            return 0.0;
        }
    });

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorInverseDistanceSummationResponseUtils::CalculateGradient(
    ModelPart &rModelPart,
    const double P,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    auto p_result = Kratos::make_shared<TensorAdaptor<double>>(rModelPart.pNodes(), Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, rModelPart.NumberOfNodes())));
    auto result_data_view = p_result->ViewData();

    IndexPartition<IndexType>(rModelPart.NumberOfNodes()).for_each([P, &result_data_view, &rModelPart, &rDistanceMatrix](const auto k) {
        double& r_value = result_data_view[k];
        r_value = 0.0;
        const double coeff = std::pow((rModelPart.NodesBegin() + k)->GetValue(SENSOR_STATUS), P - 1) * P;
        for (IndexType i = 0; i < rModelPart.NumberOfNodes(); ++i) {
            if (i != k) {
                r_value += std::pow((rModelPart.NodesBegin() + i)->GetValue(SENSOR_STATUS), P) * coeff / rDistanceMatrix.GetDistance(i, k);
            }
        }
    });

    return p_result;

    KRATOS_CATCH("");
}

} // namespace Kratos