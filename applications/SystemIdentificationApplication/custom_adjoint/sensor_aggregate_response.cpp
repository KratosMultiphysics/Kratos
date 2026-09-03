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
#include <numeric>


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
        double mDerivativeScale;

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
            .mDerivativeScale = 0.0,
            .mScaledErrors = std::vector<double>(mpImpl->mSensors.size())});
        Impl::Cache& r_cache = mpImpl->mMaybeCache.value();

        // Precompute required components.
        IndexPartition<>(mpImpl->mSensors.size()).for_each(
            [&r_cache, &rModelPart, this] (std::size_t i_sensor) -> void {
                SensorResponse& r_sensor = *mpImpl->mSensors[i_sensor];
                r_sensor.ComputeCache(rModelPart);

                const double sensor_scale = mpImpl->GetSensorScale(i_sensor);
                const double sensor_value = r_sensor.ComputeValue(rModelPart, 0);
                const double sensor_error = sensor_value - r_sensor.GetNode()->GetValue(SENSOR_MEASURED_VALUE);
                const double scaled_sensor_error = sensor_error * sensor_scale;
                const double inner_component = std::pow(
                    scaled_sensor_error,
                    mpImpl->mExponent);

                r_cache.mScaledErrors[i_sensor] = scaled_sensor_error;
                AtomicAdd(r_cache.mDerivativeScale, inner_component);
            }); // IndexPartition(mSensors.size())

        r_cache.mDerivativeScale = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(r_cache.mDerivativeScale);
        r_cache.mDerivativeScale = (1.0 / mpImpl->mExponent) * std::pow(
            r_cache.mDerivativeScale,
            (1.0 / mpImpl->mExponent) - 1.0);
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
                    return rp_left->GetContainingElement().value()->Id() < rp_right->GetContainingElement().value()->Id();
                });
        KRATOS_CATCH("");
}


double SensorAggregateResponse::ComputeValue(
    const ModelPart& rModelPart,
    int iBuffer) const {
        KRATOS_TRY
            const Impl::Cache& r_cache = mpImpl->mMaybeCache.value();
            const double inner_sum = std::accumulate(
                r_cache.mScaledErrors.begin(),
                r_cache.mScaledErrors.end(),
                0.0,
                [this] (double sum, double scaled_error) -> double {
                    return sum + std::pow(
                        scaled_error,
                        mpImpl->mExponent);
                });

            return std::pow(inner_sum, 1.0 / mpImpl->mExponent);
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


void SensorAggregateResponse::ComputeDerivative(
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
                    return rp_sensor->GetContainingElement().value()->Id() < id_element;
                });
            const auto it_sensor_end = std::upper_bound(
                it_sensor_begin,
                mpImpl->mSensors.end(),
                rElement.Id(),
                [] (Element::IndexType id_element, const SensorResponse::Pointer& rp_sensor) -> bool {
                    return id_element < rp_sensor->GetContainingElement().value()->Id();
                });

            const std::size_t i_sensor_begin = std::distance(
                mpImpl->mSensors.begin(),
                it_sensor_begin);
            const std::size_t i_sensor_end = std::distance(
                mpImpl->mSensors.begin(),
                it_sensor_end);

            // Loop over relevant sensors and sum their contributions up.
            Vector sensor_derivative_buffer;
            std::array<std::vector<IAdjoint::DynamicVariable>,2> sensor_variable_buffers; // <== state_vars, design_vars
            std::vector<std::size_t> sensor_derivative_map;
            const Impl::Cache& r_cache = mpImpl->mMaybeCache.value();

            for (const std::size_t i_sensor : std::views::iota(i_sensor_begin, i_sensor_end)) {
                const SensorResponse::Pointer& rp_sensor = mpImpl->mSensors[i_sensor];

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
                const double sensor_scale = rp_sensor->GetNode()->GetValue(SENSOR_NORMALIZATION_FACTOR);
                const double derivative_scale =
                    mpImpl->mExponent
                    * std::pow(r_cache.mScaledErrors[i_sensor], mpImpl->mExponent - 1)
                    * sensor_scale;

                // Compute, map and reduce relevant derivatives.
                for (auto& r_sensor_variables : sensor_variable_buffers) {
                    KRATOS_TRY
                        rp_sensor->ComputeDerivative(
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
            } // for i_sensor in range(i_sensor_begin, i_sensor_end)

            rOutput *= mpImpl->mMaybeCache.value().mDerivativeScale;
        KRATOS_CATCH("")
}


void SensorAggregateResponse::ComputeDerivative(
    Vector& rOutput,
    const Condition&,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int iBuffer) const {
        rOutput = ZeroVector(Variables.size());
}


} // namespace Kratos
