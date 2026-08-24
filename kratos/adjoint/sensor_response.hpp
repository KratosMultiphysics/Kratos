//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Kratos Core Includes ---
#include "adjoint/response_function.hpp"
#include "tensor_adaptors/tensor_adaptor.h"


namespace Kratos {


class KRATOS_API(KRATOS_CORE) SensorResponse : public ResponseFunction {
public:
    KRATOS_CLASS_POINTER_DEFINITION(SensorResponse);

    inline static constexpr double DefaultErrorThreshold = 1e-16;

    /**
     * @brief Construct a new SensorResponse object
     * @details         This construct a new sensor with a given name, weigh and a node. The node
     *                  position should reflect the location of the sensor. It does not have to necessarily
     *                  coincide with a mesh node. Usually these nodes are created in a separate model part
     *                  and not as a part of the mesh.
     * @param DesignVariableTypes   List of design variables.
     * @param rName                 Name of the sensor. Needs to be unique.
     * @param pNode                 Node which represents the location of the sensor and data value container.
     * @param pElement              Element in which the sensor is located.
     * @param Weight                Weight of the sensor.
     * @param ErrorThreshold        Absolute error threshold below which SENSOR_ERROR is treated as zero.
     */
    SensorResponse(
        std::span<const Variable<double>* const> DesignVariableTypes,
        const std::string& rName,
        Node::Pointer pNode,
        Element::Pointer pElement,
        const double Weight,
        const double ErrorThreshold = DefaultErrorThreshold);

    /**
     * @brief Creates a new sensor attached to a domain model part given.
     * @details This static method is used to construct the sensors on a domain model part
     *          with settings passed in using SensorParameters.
     *
     * @param rDomainModelPart      Domain model part, from which sensor will compute its values, sensitivities, etc.
     * @param rSensorModelPart      Model part, which will be used to create the node for the new sensor. This node will be having the location of the sensor.
     * @param Id                    Id of the created sensor.
     * @param SensorParameters      Parameters required to construct the sensor.
     * @return SensorResponse::Pointer      A new sensor.
     */
    virtual SensorResponse::Pointer Create(
        const ModelPart& rDomainModelPart,
        const ModelPart& rSensorModelPart,
        IndexType Id,
        Parameters SensorParameters) const {
            KRATOS_ERROR << "Trying to create a default \"SensorResponse\" which is not allowed. Please create specific sensor." << std::endl;
            return SensorResponse::Pointer();
    }

    virtual Parameters GetDefaultParameters() const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    /**
     * @brief Get the SensorResponse parameters.
     *
     * This method is used retrieve sensor parameters used in this sensor.
     * These parameters are used to reconstruct the sensor from scratch
     * if required.
     *
     * @return const Parameters
     */
    virtual Parameters GetSensorParameters() const;

    /**
     * @brief Get the Type of the sensor
     *
     * @return std::string
     */
    const std::string& GetName() const noexcept {
        return mName;
    }

    /// @brief Get the sensor location.
    intrusive_ptr<const Node> GetNode() const noexcept {
        return mpNode;
    }

    /// @brief Get the sensor location.
    Node::Pointer GetNode() noexcept {
        return mpNode;
    }

    /**
     * @brief Get the Weight of the sensor.
     *
     * @return double       The weight of the sensor .
     */
    double GetWeight() const noexcept {
        return mWeight;
    }

    /**
     * @brief Get the Error Threshold of the sensor.
     *
     * @return double       The absolute error threshold below which SENSOR_ERROR is treated as zero.
     */
    double GetErrorThreshold() const noexcept {
        return mErrorThreshold;
    }

    /**
     * @brief Adds given tensor adaptor
     *
     * @throws If an tensor adaptor is found already under the given name.
     *
     * @param rTensorAdaptorName       Name of the tensor adaptor.
     * @param pTensorAdaptor           Tensor adaptor pointer to be added.
     */
    void AddTensorAdaptor(
        const std::string& rTensorAdaptorName,
        TensorAdaptor<double>::Pointer pTensorAdaptor);

    /**
     * @brief Get the tensor adaptor for specified tensor adaptor name
     *
     * @throws If the @ref rTensorAdaptorName is not found in the map of tensor adaptors.
     *
     * @param rTensorAdaptorName                Tensor adaptor name
     * @return TensorAdaptor<double>::Pointer   Tensor adaptor
     */
    TensorAdaptor<double>::Pointer GetTensorAdaptor(const std::string& rTensorAdaptorName) const;

    /// @brief Clear tensor adaptors from the sensors.
    void ClearTensorAdaptors();

    virtual std::string Info() const;

    virtual void PrintInfo(std::ostream& rOStream) const;

    virtual void PrintData(std::ostream& rOStream) const;

protected:
    const Element& GetContainingElement() const {
        return *mpElement;
    }

    Element& GetContainingElement() {
        return *mpElement;
    }

private:
    std::string mName;

    Node::Pointer mpNode;

    Element::Pointer mpElement;

    double mWeight;

    double mErrorThreshold;

    std::unordered_map<
        std::string,
        TensorAdaptor<double>::Pointer
    > mTensorAdaptorsMap;
}; // class SensorResponse


} // namespace Kratos
