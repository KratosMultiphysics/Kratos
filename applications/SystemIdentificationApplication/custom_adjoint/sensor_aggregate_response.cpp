//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Kratos SysId Includes ---
#include "custom_adjoint/sensor_aggregate_response.hpp"
#include "system_identification_application_variables.h"

// --- Kratos Core Includes ---
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

// --- STL Includes ---
#include <optional>
#include <ranges>


namespace Kratos {


struct SensorAggregateResponse::Impl {
    std::vector<SensorResponse::Pointer> mSensors;

    double mExponent;

    /// @brief Components required for computing the value and derivative of this response.
    /// @details These components are computed in @ref InitializeSolutionStep and erased in
    ///          @ref FinalizeSolutionStep. They rely on specific variables being stored in
    ///          the nodes belonging to the member sensors.
    ///
    ///          Notation in the docs of the member variables:
    ///          - @f$v_i@f$: @ref SensorResponse::ComputeValue
    ///          - @f$m_i@f$: @p SensorResponse::GetNode::GetValue(SENSOR_MEASURED_VALUE)
    ///          - @f$s_i@f$: scaling factor for the given sensor @p SensorResponse::GetNode::GetValue(SENSOR_NORMALIZATION_FACTOR)
    struct Cache {
        /// @f[
        ///     \sum_i \left( (v_i - m_i) s_i \right)^p
        /// @f]
        double mValue;

        /// @details Stores the normalized and weighted error for each sensor.
        ///          @f[
        ///              (v_i - m_i) s_i
        ///          @f]
        std::vector<double> mScaledErrors;
    };
    std::optional<Cache> mMaybeCache;

    [[nodiscard]] double GetSensorScale(std::size_t iSensor) const {
        return mSensors[iSensor]->GetNode()->GetValue(SENSOR_NORMALIZATION_FACTOR);
    }
}; // struct SensorAggregateResponse::Impl


SensorAggregateResponse::SensorAggregateResponse() noexcept
    : SensorAggregateResponse(2ul)
{}


SensorAggregateResponse::SensorAggregateResponse(std::size_t Exponent)
    : mpImpl(new Impl {
        .mSensors = {},
        .mExponent = static_cast<double>(Exponent),
        .mMaybeCache = {}}) {
            KRATOS_ERROR_IF(Exponent <= 0ul || Exponent % 2ul)
                << std::format("exponent ({}) must be a positive even integer", Exponent);
}


SensorAggregateResponse::SensorAggregateResponse(SensorAggregateResponse&&) noexcept = default;


SensorAggregateResponse::SensorAggregateResponse(const SensorAggregateResponse& rRhs)
    : SensorAggregateResponse() {
        *mpImpl = *rRhs.mpImpl;
}


SensorAggregateResponse& SensorAggregateResponse::operator=(SensorAggregateResponse&& rRhs) noexcept = default;


SensorAggregateResponse& SensorAggregateResponse::operator=(const SensorAggregateResponse& rRhs) {
    *mpImpl = *rRhs.mpImpl;
    return *this;
}


SensorAggregateResponse::~SensorAggregateResponse() = default;


void SensorAggregateResponse::ComputeCache(const ModelPart& rModelPart) {
    KRATOS_TRY
        // Construct a new cache.
        mpImpl->mMaybeCache.emplace(Impl::Cache {
            .mValue = 0.0,
            .mScaledErrors = std::vector<double>(mpImpl->mSensors.size())});
        Impl::Cache& r_cache = mpImpl->mMaybeCache.value();

        // Precompute required components.
        IndexPartition<>(mpImpl->mSensors.size()).for_each(
            [&r_cache, &rModelPart, this] (std::size_t i_sensor) -> void {
                SensorResponse& r_sensor = *mpImpl->mSensors[i_sensor];
                r_sensor.ComputeCache(rModelPart);

                const double sensor_scale = mpImpl->GetSensorScale(i_sensor);
                const double sensor_value = r_sensor.ComputeInnerValue(
                    r_sensor.GetContainingElement(),
                    rModelPart.GetProcessInfo(),
                    0);
                const double sensor_error = sensor_value - r_sensor.GetNode()->GetValue(SENSOR_MEASURED_VALUE);
                const double scaled_sensor_error = sensor_error * sensor_scale;
                r_cache.mScaledErrors[i_sensor] = scaled_sensor_error;

                const double inner_component = std::pow(
                    scaled_sensor_error,
                    mpImpl->mExponent);
                AtomicAdd(r_cache.mValue, inner_component);
            }); // IndexPartition(mSensors.size())

        r_cache.mValue = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_cache.mValue);
    KRATOS_CATCH("")
}


void SensorAggregateResponse::ClearCache() {
    KRATOS_TRY
        block_for_each(
            mpImpl->mSensors,
            [] (SensorResponse::Pointer& rp_sensor) -> void {
                rp_sensor->ClearCache();
            });
    KRATOS_CATCH("")
    mpImpl->mMaybeCache.reset();
}


void SensorAggregateResponse::AddSensors(std::span<const SensorResponse::Pointer> Sensors) {
        KRATOS_TRY
            // Check requirements and assign design variable types.
            block_for_each(
                Sensors,
                [this] (const auto& rp_sensor) -> void {
                    KRATOS_ERROR_IF_NOT(rp_sensor->GetNode()->Has(SENSOR_MEASURED_VALUE))
                        << rp_sensor->Name() << " requires SENSOR_MEASURED_VALUE";
                    KRATOS_ERROR_IF_NOT(rp_sensor->GetNode()->Has(SENSOR_NORMALIZATION_FACTOR))
                        << rp_sensor->Name() << " requires SENSOR_NORMALIZATION_FACTOR";
                    rp_sensor->SetDesignVariableTypes(this->GetDesignVariableTypes());
                });

            // Insert all sensors.
            mpImpl->mSensors.insert(
                mpImpl->mSensors.end(),
                Sensors.begin(),
                Sensors.end());

            // Sort sensors by the IDs of the elements they are located in.
            std::sort(
                mpImpl->mSensors.begin(),
                mpImpl->mSensors.end(),
                [] (const auto& rp_left, const auto& rp_right) -> bool {
                    return rp_left->GetContainingElement().Id() < rp_right->GetContainingElement().Id();
                });
        KRATOS_CATCH("");
}


double SensorAggregateResponse::ComputeInnerValue(
    const Element& rElement,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_TRY
            // Find which sensors are located in the provided element.
            const auto it_sensor_begin = std::lower_bound(
                mpImpl->mSensors.begin(),
                mpImpl->mSensors.end(),
                rElement.Id(),
                [] (const SensorResponse::Pointer& rp_sensor, Element::IndexType id_element) -> bool {
                    return rp_sensor->GetContainingElement().Id() < id_element;
                });
            const auto it_sensor_end = std::upper_bound(
                it_sensor_begin,
                mpImpl->mSensors.end(),
                rElement.Id(),
                [] (Element::IndexType id_element, const SensorResponse::Pointer& rp_sensor) -> bool {
                    return id_element < rp_sensor->GetContainingElement().Id();
                });
            const auto sensor_range = std::ranges::subrange(
                it_sensor_begin,
                it_sensor_end);

            // Loop over relevant sensors and sum their contributions up.
            const double element_contribution = std::accumulate(
                sensor_range.begin(),
                sensor_range.end(),
                0.0,
                [rElement, rProcessInfo, iBuffer, this] (double sum, const SensorResponse::Pointer& rp_sensor) -> double {
                    const double virtual_sensor_value = rp_sensor->ComputeInnerValue(
                        rElement,
                        rProcessInfo,
                        iBuffer);
                    const double measured_sensor_value = rp_sensor->GetNode()->GetValue(SENSOR_MEASURED_VALUE);
                    const double scale = rp_sensor->GetNode()->GetValue(SENSOR_NORMALIZATION_FACTOR);
                    const double scaled_error = (virtual_sensor_value - measured_sensor_value) * scale;
                    return sum + std::pow(scaled_error, mpImpl->mExponent);
                });

            return element_contribution;
        KRATOS_CATCH("");
}


void SensorAggregateResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Element& rElement,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rVariables.clear();
            std::vector<IAdjoint::DynamicVariable> sensor_variables;
            for (const SensorResponse::Pointer& rp_sensor : mpImpl->mSensors) {
                sensor_variables.clear();
                rp_sensor->GetStateVariables(sensor_variables, rElement, rProcessInfo);
                rVariables.insert(rVariables.end(), sensor_variables.begin(), sensor_variables.end());
            }

            std::sort(
                rVariables.begin(),
                rVariables.end());
            rVariables.erase(
                std::unique(
                    rVariables.begin(),
                    rVariables.end()),
                rVariables.end());
        KRATOS_CATCH("")
}


void SensorAggregateResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition&,
    const ProcessInfo&) const {
        rVariables.clear();
}


void SensorAggregateResponse::ComputeInnerDerivative(
    Vector& rOutput,
    const Element& rElement,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";

        KRATOS_TRY
            // Wipe the output.
            rOutput = ZeroVector(Variables.size());

            // Find which sensors are located in the provided element.
            const auto it_sensor_begin = std::lower_bound(
                mpImpl->mSensors.begin(),
                mpImpl->mSensors.end(),
                rElement.Id(),
                [] (const SensorResponse::Pointer& rp_sensor, Element::IndexType id_element) -> bool {
                    return rp_sensor->GetContainingElement().Id() < id_element;
                });
            const auto it_sensor_end = std::upper_bound(
                it_sensor_begin,
                mpImpl->mSensors.end(),
                rElement.Id(),
                [] (Element::IndexType id_element, const SensorResponse::Pointer& rp_sensor) -> bool {
                    return id_element < rp_sensor->GetContainingElement().Id();
                });
            const auto sensor_range = std::ranges::subrange(
                it_sensor_begin,
                it_sensor_end);

            // Loop over relevant sensors and sum their contributions up.
            Vector sensor_derivative_buffer;
            std::array<std::vector<IAdjoint::DynamicVariable>,2> sensor_variable_buffers; // <== state_vars, design_vars
            std::vector<std::size_t> sensor_derivative_map;

            for (const SensorResponse::Pointer& rp_sensor : sensor_range) {
                // Collect the list of relevant variables for the current sensor.
                rp_sensor->GetStateVariables(
                    sensor_variable_buffers.front(),
                    rElement,
                    rProcessInfo);
                rp_sensor->GetDesignVariables(
                    sensor_variable_buffers.back(),
                    rElement,
                    rProcessInfo);
                for (auto& r_variable_buffer : sensor_variable_buffers)
                    r_variable_buffer.erase(
                        std::remove_if(
                            r_variable_buffer.begin(),
                            r_variable_buffer.end(),
                            [&Variables] (const IAdjoint::DynamicVariable& r_variable) -> bool {
                                return std::find(
                                    Variables.begin(),
                                    Variables.end(),
                                    r_variable) == Variables.end();}),
                        r_variable_buffer.end());

                if (sensor_variable_buffers.front().empty() && sensor_variable_buffers.back().empty())
                    continue;

                // Compute the sensor's output value.
                const double virtual_sensor_value = rp_sensor->ComputeInnerValue(
                    rElement,
                    rProcessInfo,
                    iBuffer);
                const double measured_value = rp_sensor->GetNode()->GetValue(SENSOR_MEASURED_VALUE);
                const double sensor_scale = rp_sensor->GetNode()->GetValue(SENSOR_NORMALIZATION_FACTOR);
                const double derivative_scale =
                    mpImpl->mExponent
                    * std::pow(virtual_sensor_value - measured_value, mpImpl->mExponent - 1)
                    * std::pow(sensor_scale, mpImpl->mExponent);

                // Compute, map and reduce relevant derivatives.
                for (auto& r_sensor_variables : sensor_variable_buffers) {
                    KRATOS_TRY
                        rp_sensor->ComputeInnerDerivative(
                            sensor_derivative_buffer,
                            rElement,
                            sensor_variable_buffers.front(),
                            rProcessInfo,
                            iBuffer);
                    KRATOS_CATCH("during the derivative evaluation of sensor " + rp_sensor->Name());
                    for (std::size_t i_variable=0ul; i_variable<sensor_variable_buffers.front().size(); ++i_variable) {
                        const std::size_t i_output = std::distance(
                            Variables.begin(),
                            std::find(
                                Variables.begin(),
                                Variables.end(),
                                r_sensor_variables[i_variable]));
                        rOutput[i_output] += derivative_scale * sensor_derivative_buffer[i_variable];
                    } // for i_variable
                } // for r_sensor_variables in sensor_variable_buffers
            } // for rp_sensor in sensor_range
        KRATOS_CATCH("")
}


void SensorAggregateResponse::ComputeInnerDerivative(
    Vector& rOutput,
    const Condition&,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int iBuffer) const {
        KRATOS_ERROR_IF_NOT(iBuffer == 0)
            << "requested buffer " << iBuffer << ", which is not supported";

        KRATOS_TRY
            // Sensors aren't supposed to provide anything on conditions.
            for (const IAdjoint::DynamicVariable& r_variable : Variables) {
                // Check whether the requested variable is a design variable.
                const auto& r_design_variable_types = this->GetDesignVariableTypes();
                const auto it_design_variable_type = std::find_if(
                    r_design_variable_types.begin(),
                    r_design_variable_types.end(),
                    [&r_variable] (const Variable<double>* p_variable_type) -> bool {
                        return p_variable_type->Key() == r_variable;
                    });
                const bool is_design_variable = it_design_variable_type != r_design_variable_types.end();

                // Error if the requested variable is neither a state nor a design variable.
                KRATOS_ERROR_IF_NOT(is_design_variable) << "unsupported variable " << r_variable.Name();
            }

            rOutput = ZeroVector(Variables.size());
        KRATOS_CATCH("")
}


} // namespace Kratos
