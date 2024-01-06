//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

// System includes
#include <algorithm>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "sensor_distance_matrix.h"

namespace Kratos {

SensorDistanceMatrix::SensorDistanceMatrix(const std::vector<Sensor::Pointer>& rSensors)
    : mSensors(rSensors)
{
    KRATOS_TRY

    const auto n = rSensors.size();

    if (n > 0) {
        mCompressedDistanceMatrix.resize(static_cast<IndexType>(n * (n - 1) / 2));
        IndexPartition<IndexType>(mCompressedDistanceMatrix.size()).for_each([&](const auto Index) {
            for (IndexType i = 0; i < n; ++i) {
                const auto j = Index + static_cast<IndexType>((i + 2) * (i + 1) / 2) - n * i;
                if (i < j && j >= 0 && j < n) {
                    mCompressedDistanceMatrix[Index] = norm_2(rSensors[i]->GetLocation() - rSensors[j]->GetLocation());
                    break;
                }
            }
        });
    } else {
        mCompressedDistanceMatrix.clear();
    }

    KRATOS_CATCH("");
}

double SensorDistanceMatrix::GetDistance(
    const Sensor& rSensor1,
    const Sensor& rSensor2) const
{
    if (&rSensor1 == &rSensor2) {
        return 0.0;
    } else {
        auto p_itr_1 = std::find_if(mSensors.begin(), mSensors.end(), [&rSensor1](const auto& rSensor) { return &rSensor1 == rSensor.get(); });
        auto p_itr_2 = std::find_if(mSensors.begin(), mSensors.end(), [&rSensor2](const auto& rSensor) { return &rSensor2 == rSensor.get(); });

        if (p_itr_1 != mSensors.end() && p_itr_2 != mSensors.end()) {
            const auto i = std::distance(mSensors.begin(), p_itr_1);
            const auto j = std::distance(mSensors.begin(), p_itr_2);

            IndexType local_i, local_j;
            if (i < j) {
                local_i = i;
                local_j = j;
            } else {
                local_i = j;
                local_j = i;
            }

            return mCompressedDistanceMatrix[mSensors.size() * local_i + local_j - static_cast<IndexType>((local_i + 2) * (local_i + 1) / 2)];
        } else {
            return -1.0;
        }
    }
}

Matrix SensorDistanceMatrix::GetDistance(
    const std::vector<Sensor::Pointer>& rSensorsList1,
    const std::vector<Sensor::Pointer>& rSensorsList2) const
{
    KRATOS_TRY

    Matrix result(rSensorsList1.size(), rSensorsList2.size());

    IndexPartition<IndexType>(rSensorsList1.size() * rSensorsList2.size()).for_each([&](const auto Index) {
        const auto& r_sensor_1 = *rSensorsList1[static_cast<IndexType>(Index / rSensorsList2.size())];
        const auto& r_sensor_2 = *rSensorsList2[Index % rSensorsList2.size()];
        result.data().begin()[Index] = this->GetDistance(r_sensor_1, r_sensor_2);
    });

    return result;

    KRATOS_CATCH("");
}

Matrix SensorDistanceMatrix::MaxDistances(
    const Matrix& rM1,
    const Matrix& rM2)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rM1.size1() == rM2.size1() && rM1.size2() == rM2.size2())
        << "Matrix shape mismatch [ rM1.shape() = (" << rM1.size1() << ", " << rM1.size2()
        << "), rM2.shape() = (" << rM2.size1() << ", " << rM2.size2() << ") ].\n";

    Matrix result(rM1.size1(), rM1.size2());
    IndexPartition<IndexType>(rM1.size1() * rM1.size2()).for_each([&result, &rM1, &rM2](const auto Index){
        result.data().begin()[Index] = std::max(rM1.data().begin()[Index], rM2.data().begin()[Index]);
    });

    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/