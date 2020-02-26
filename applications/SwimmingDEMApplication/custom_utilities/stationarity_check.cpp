//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas
//

// System includes

// External includes

// Project includes
#include "stationarity_check.h"

namespace Kratos {

bool FlowStationarityCheck::AssessStationarity()
{
    // Loop over all nodes:
    ModelPart::NodesContainerType& nodes_array = mrModelPart.Nodes();
    const int num_elem = static_cast<int>(nodes_array.size());

    // Compute the error estimate per element
    double average_pressure_change_rate = 0.0;
    const double time_step_reciprocal = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];
    //TODO: add error if delta_time is 0

    #pragma omp parallel for reduction(+:average_pressure_change_rate)
    for (int i = 0; i < num_elem; ++i){
        auto it_node = nodes_array.begin() + i;
        const double old_pressure = it_node->FastGetSolutionStepValue(PRESSURE, 1);
        const double pressure = it_node->FastGetSolutionStepValue(PRESSURE);
        average_pressure_change_rate += time_step_reciprocal * std::abs(pressure - old_pressure);
    }

    mCurrentPressureRate = average_pressure_change_rate;

    if (this->GetTransienceMeasure() < mTolerance){
        mCharacteristicPressureRate = (mAveragingStep * mCharacteristicPressureRate + average_pressure_change_rate) / (mAveragingStep + 1);
        ++mAveragingStep;
        return true;
    }

    else {
        mCharacteristicPressureRate = (mAveragingStep * mCharacteristicPressureRate + average_pressure_change_rate) / (mAveragingStep + 1);
        ++mAveragingStep;
        return false;
    }
}

double FlowStationarityCheck::GetCharacteristicPressureDerivative()
{
    return mCharacteristicPressureRate;
}

double FlowStationarityCheck::GetCurrentPressureDerivative()
{
    return mCurrentPressureRate;
}

double FlowStationarityCheck::GetTolerance()
{
    return mTolerance;
}

double FlowStationarityCheck::GetTransienceMeasure()
{
    return mCurrentPressureRate / mCharacteristicPressureRate;
}

}  // namespace Kratos.


