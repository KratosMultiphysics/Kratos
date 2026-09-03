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

// --- STL Includes ---
#include <optional>


namespace Kratos {


/// @brief Response function defined over a single physical point.
/// @details Sensors must be located inside (or on the boundary of) a single element. Therefore,
///          sensors are expected to only provide contributions for @ref Element "elements",
///          but not @ref Condition "conditions".
/// @ingroup adjoints
/// @todo Relax the requirement of being tied to a single point (i.e.: support sensors over a geometry). @matekelemen
class KRATOS_API(KRATOS_CORE) SensorResponse : public ResponseFunction {
public:
    KRATOS_CLASS_POINTER_DEFINITION(SensorResponse);

    inline static constexpr double DefaultErrorThreshold = 1e-16;

    SensorResponse() noexcept = default;

    /**
     * @brief Construct a new SensorResponse object
     * @details         This construct a new sensor with a given name, weigh and a node. The node
     *                  position should reflect the location of the sensor. It does not have to necessarily
     *                  coincide with a mesh node. Usually these nodes are created in a separate model part
     *                  and not as a part of the mesh.
     * @param DesignVariableTypes   List of design variables.
     * @param rName                 Name of the sensor. Needs to be unique.
     * @param pNode                 Node which represents the location of the sensor and data value container.
     * @param pMaybeElement         Element in which the sensor is located (if the element is on this rank).
     */
    SensorResponse(
        std::span<const Variable<double>* const> DesignVariableTypes,
        const std::string& rName,
        Node::Pointer pNode,
        std::optional<intrusive_ptr<const Element>> pMaybeElement);

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
    [[nodiscard]] virtual SensorResponse::Pointer Create(
        [[maybe_unused]] const ModelPart& rDomainModelPart,
        [[maybe_unused]] const ModelPart& rSensorModelPart,
        [[maybe_unused]] IndexType Id,
        [[maybe_unused]] Parameters SensorParameters) const {
            KRATOS_ERROR << "Trying to create a default \"SensorResponse\" which is not allowed. Please create specific sensor." << std::endl;
            return SensorResponse::Pointer();
    }

    /// @name Adjoint Interface
    /// @{

    /// @copydoc ResponseFunction::GetStateVariables(std::vector<IAdjoint::DynamicVariable>&,const Condition&,const ProcessInfo&) const
    void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const final override;

    using ResponseFunction::GetStateVariables;

    /// @copydoc ResponseFunction::GetDesignVariables(std::vector<IAdjoint::DynamicVariable>&,const Condition&,const ProcessInfo&) const
    void GetDesignVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const final override;

    using ResponseFunction::GetDesignVariables;

    /// @copydoc ResponseFunction::ComputeDerivative(Vector&,const Condition&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeDerivative(
        Vector& rOutput,
        const Condition& rCondition,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const final override;

    using ResponseFunction::ComputeDerivative;

    /// @}

    virtual Parameters GetDefaultParameters() const {
        KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
    }

    [[nodiscard]] const std::string& Name() const noexcept {
        return mName;
    }

    /// @brief Get the sensor location.
    [[nodiscard]] intrusive_ptr<const Node> GetNode() const noexcept {
        return mpNode;
    }

    /// @brief Get the sensor location.
    [[nodiscard]] Node::Pointer GetNode() noexcept {
        return mpNode;
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
    [[nodiscard]] TensorAdaptor<double>::Pointer GetTensorAdaptor(const std::string& rTensorAdaptorName) const;

    /// @brief Clear tensor adaptors from the sensors.
    void ClearTensorAdaptors();

    [[nodiscard]] virtual std::string Info() const;

    virtual void PrintInfo(std::ostream& rOStream) const;

    virtual void PrintData(std::ostream& rOStream) const;

    [[nodiscard]] std::optional<const Element*> GetContainingElement() const {
        if (!mpMaybeElement.has_value()) return {}; else return mpMaybeElement.value().get();
    }

private:
    std::string mName;

    Node::Pointer mpNode;

    std::optional<intrusive_ptr<const Element>> mpMaybeElement;

    std::unordered_map<
        std::string,
        TensorAdaptor<double>::Pointer
    > mTensorAdaptorsMap;
}; // class SensorResponse


} // namespace Kratos
