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

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_sensors/sensor.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorDistanceMatrix
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SensorDistanceMatrix);

    ///@}
    ///@name Life cycle
    ///@{

    SensorDistanceMatrix(const std::vector<Sensor::Pointer>& rSensors);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Get the Compressed Distance Matrix
     *
     * @return const std::vector<double>&
     */
    const std::vector<double>& GetCompressedDistanceMatrix() const { return mCompressedDistanceMatrix; }

    /**
     * @brief Get the Distance between two sensors.
     *
     * This method returns distance between two sensors if they are found in the Sensors list
     * which is used to construct this SensorDistaneMatrix. If not found, then -1.0 is returned.
     *
     * @param rSensor1          Sensor 1.
     * @param rSensor2          Sensor 2.
     * @return double           Distance between two sensors if found, otherwise -1.0.
     */
    double GetDistance(
        const Sensor& rSensor1,
        const Sensor& rSensor2) const;

    /**
     * @brief Get the Distance matrix between two sensors lists.
     *
     * This returns a matrix with size (n1, n2) where n1 corresponds to number
     * of sensors in rSensorsList1 and n2 corresponds to rSensorsList2. The matrix
     * will contain -1.0 values for each entry if either component of the sensors lists
     * are not found in the sensors list which is used at the construction of this
     * SensorDistanceMatrix.
     *
     * @param rSensorsList1     Sensors list 1.
     * @param rSensorsList2     Sensors list 2.
     * @return Matrix           Matrix with distances.
     */
    Matrix GetDistance(
        const std::vector<Sensor::Pointer>& rSensorsList1,
        const std::vector<Sensor::Pointer>& rSensorsList2) const;

    /**
     * @brief Get the max distances from two matrices
     *
     * @param rM1               Matrix 1.
     * @param rM2               Matrix 2.
     * @return Matrix           Resulting matrix with maximums.
     */
    static Matrix MaxDistances(
        const Matrix& rM1,
        const Matrix& rM2);

    ///@}
private:
    ///@name Private member variables
    ///@{

    const std::vector<Sensor::Pointer> mSensors;

    std::vector<double> mCompressedDistanceMatrix;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/