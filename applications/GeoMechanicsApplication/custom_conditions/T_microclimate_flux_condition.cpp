// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    John van Esch
//                   Mohamed Nabi
//

// Application includes
#include "custom_conditions/T_microclimate_flux_condition.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "custom_utilities/ublas_utils.h"

namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TMicroClimateFluxCondition<TDim, TNumNodes>::Create(
    IndexType NewId,
    const NodesArrayType& rNodes,
    Properties::Pointer pProperties) const
{
    return make_intrusive<TMicroClimateFluxCondition>(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    if (!mIsInitialised)
    {
        this->InitializeProperties();
        mIsInitialised = true;
    }
    this->CalculateRoughness(rCurrentProcessInfo);
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(
    Vector& rRightHandSideVector,
    const Vector& rNodalTemperatures)
{
    auto temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{outer_prod(mVariables.Np, mVariables.Np) * mVariables.IntegrationCoefficient};
    auto temporary_vector = array_1d<double,TNumNodes>{prod(temporary_matrix, mVariables.rightHandSideFlux)};
    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, temporary_vector);

    auto flux_matrix = Matrix{TNumNodes, TNumNodes, 0.0};
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        flux_matrix(i, i) = mVariables.leftHandSideFlux[i];
    }
    temporary_matrix = prod(temporary_matrix, flux_matrix);
    temporary_vector = -prod(temporary_matrix, rNodalTemperatures);
    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, temporary_vector);
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddLHS(
    Matrix& rLeftHandSideMatrix)
{
    KRATOS_TRY

    auto temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{outer_prod(mVariables.Np, mVariables.Np) * mVariables.IntegrationCoefficient};

    auto flux_matrix = Matrix{TNumNodes, TNumNodes, 0.0};
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        flux_matrix(i, i) = mVariables.leftHandSideFlux[i];
    }

    temporary_matrix = prod(temporary_matrix, flux_matrix);

    GeoElementUtilities::
        AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, temporary_matrix);

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim,TNumNodes>::CalculateIntegrationCoefficient(
    const Matrix& Jacobian,
    double Weight)
{
    if (TDim == 2)
    {
        return MathUtils<>::Norm(UBlasUtils::MakeVector({Jacobian(0, 0), Jacobian(1, 0)})) * Weight;
    }
    else if (TDim == 3)
    {
        const auto vec1 = UBlasUtils::MakeVector({Jacobian(0, 0), Jacobian(1, 0), Jacobian(2, 0)});
        const auto vec2 = UBlasUtils::MakeVector({Jacobian(0, 1), Jacobian(1, 1), Jacobian(2, 1)});
        return MathUtils<>::Norm(MathUtils<>::CrossProduct(vec1, vec2)) * Weight;
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRoughness(
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto time_step_size = rCurrentProcessInfo.GetValue(DELTA_TIME);

    const auto& r_geom = this->GetGeometry();
    const auto current_air_temperature = r_geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto current_wind_speed = std::max(1.0e-3, r_geom[0].FastGetSolutionStepValue(WIND_SPEED));

    constexpr auto roughness_layer_height = 10.0;
    constexpr auto roughness_layer_resistance = 30.0;
    constexpr auto von_neuman_coefficient = 0.4;
    constexpr auto measurement_height = 10.0;
    constexpr auto roughness_height = 1.0;
    constexpr auto gravitational_acceleration = 9.81;

    const auto previous_roughness_temperature = mRoughnessTemperature;
    mRoughnessTemperature = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto initial_soil_temperature = r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

        auto surface_roughness_factor = 0.0;

        // Eq 5.29
        const auto richardson_bulk_modulus = 2.0 * gravitational_acceleration * measurement_height / (current_air_temperature +
            previous_roughness_temperature + 546.3) * (current_air_temperature - previous_roughness_temperature) / (current_wind_speed * current_wind_speed);

        // Eq 5.25
        const auto friction_drag_coefficient = von_neuman_coefficient / std::log(measurement_height / roughness_height);

        auto cof = 0.0;
        if (previous_roughness_temperature >= current_air_temperature) {
            // Eq 5.27
            cof = richardson_bulk_modulus / (1.0 + 75.0 * friction_drag_coefficient * friction_drag_coefficient *
                std::sqrt(measurement_height / roughness_height * std::abs(richardson_bulk_modulus)));
            surface_roughness_factor = 1.0 - 15.0 * cof;
        }
        else {
            // Eq 5.28
            cof = std::sqrt(1.0 + 5.0 * richardson_bulk_modulus);
            surface_roughness_factor = 1.0 / (1.0 + 15.0 * richardson_bulk_modulus * cof);
        }

        const auto c = roughness_layer_resistance * roughness_layer_height + time_step_size + time_step_size * current_wind_speed * roughness_layer_resistance *
            surface_roughness_factor * friction_drag_coefficient * friction_drag_coefficient;
        const auto current_roughness_temperature = (roughness_layer_resistance * roughness_layer_height * previous_roughness_temperature + time_step_size *
            initial_soil_temperature + time_step_size * current_wind_speed * roughness_layer_resistance * surface_roughness_factor *
            friction_drag_coefficient * friction_drag_coefficient * current_air_temperature) / c;
        mRoughnessTemperature += current_roughness_temperature / TNumNodes;
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNodalFluxes(
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = this->GetGeometry();
    const auto time_step_size = rCurrentProcessInfo.GetValue(DELTA_TIME);

    constexpr auto air_density = 1.18;
    constexpr auto air_heat_capacity = 1004.67;
    constexpr auto roughness_layer_resistance = 30.0;
    constexpr auto latent_evaporation_heat = 2.45e6;
    constexpr auto water_density = 1e3;
    constexpr auto psychometric_constant = 0.63;
    constexpr auto surface_resistance = 30.0;

    const auto previous_storage = mWaterStorage;
    const auto previous_radiation = mNetRadiation;
    mWaterStorage = 0.0;
    mNetRadiation = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto atmospheric_temperature = r_geom[i].FastGetSolutionStepValue(AIR_TEMPERATURE);
        const auto incoming_radiation = r_geom[i].FastGetSolutionStepValue(SOLAR_RADIATION);
        const auto humidity = r_geom[i].FastGetSolutionStepValue(AIR_HUMIDITY);
        const auto precipitation = r_geom[i].FastGetSolutionStepValue(PRECIPITATION);
        const auto wind_speed = r_geom[i].FastGetSolutionStepValue(WIND_SPEED);

        // Eq 5.22
        const auto sensible_heat_flux_left = air_heat_capacity * air_density / roughness_layer_resistance;
        const auto sensible_heat_flux_right = -air_heat_capacity * air_density * mRoughnessTemperature / roughness_layer_resistance;

        // Eq 5.35
        const auto atmospheric_resistance = 1.0 / (0.007 + 0.0056 * wind_speed);

        // Eq. 5.12
        const auto saturated_vapor_pressure = 6.11 * std::exp(17.27 * atmospheric_temperature / (atmospheric_temperature + 237.3));

        // Eq 5.13
        const auto vapor_pressure_increment = 4098.0 * saturated_vapor_pressure / (std::pow((atmospheric_temperature + 237.3), 2.0));

        // Eq 5.14
        const auto actual_vapor_pressure = humidity / 100.0 * saturated_vapor_pressure;

        const auto initial_soil_temperature = r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);
        const auto net_radiation = CalculateNetRadiation(incoming_radiation, atmospheric_temperature, initial_soil_temperature);

        // Eq 5.20
        const auto surface_heat_storage = mFirstCoverStorageCoefficient * net_radiation + mSecondCoverStorageCoefficient * (net_radiation - previous_radiation) /
            time_step_size + mThirdCoverStorageCoefficient;

        // Eq 5.34
        auto latent_heat_flux = (vapor_pressure_increment * (net_radiation + mBuildEnvironmentRadiation - surface_heat_storage) + air_heat_capacity * air_density *
            (saturated_vapor_pressure - actual_vapor_pressure) / atmospheric_resistance) / (vapor_pressure_increment + psychometric_constant *  // division is replaced to * based on (3.34)
                (1.0 + surface_resistance / atmospheric_resistance));   //Where is G (5.34)?
        latent_heat_flux = std::max(latent_heat_flux, 0.0);

        auto actual_evaporation = 0.0;
        auto actual_precipitation = 0.0;

        // Eq 5.36
        const auto potential_evaporation = latent_heat_flux / (water_density * latent_evaporation_heat);
        const auto potential_storage = previous_storage + time_step_size * (precipitation - potential_evaporation);
        if (potential_storage > mMaximalStorage)
        {
            actual_evaporation = potential_evaporation;
            actual_precipitation = (mMaximalStorage - previous_storage) / time_step_size + actual_evaporation;
        }
        else if (potential_storage < mMinimalStorage)
        {
            actual_precipitation = precipitation;
            actual_evaporation = (previous_storage - mMinimalStorage) / time_step_size + actual_precipitation;
        }
        else
        {
            actual_evaporation = potential_evaporation;
            actual_precipitation = precipitation;
        }
        const auto actual_storage = previous_storage + time_step_size * (actual_precipitation - actual_evaporation);
        latent_heat_flux = actual_evaporation * water_density * latent_evaporation_heat;

        // Eq 5.31
        const auto subsurface_heat_flux = net_radiation - sensible_heat_flux_right - latent_heat_flux + mBuildEnvironmentRadiation - surface_heat_storage;

        mNetRadiation += net_radiation / TNumNodes;
        mWaterStorage += actual_storage / TNumNodes;
        mVariables.leftHandSideFlux[i] = sensible_heat_flux_left;
        mVariables.rightHandSideFlux[i] = subsurface_heat_flux;
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLocalSystem(
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rLeftHandSideMatrix = Matrix{TNumNodes, TNumNodes, 0.0};
    rRightHandSideVector = Vector{TNumNodes, 0.0};

    // Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        JContainer[i].resize(TDim, LocalDim, false);
    Geom.Jacobian(JContainer, this->GetIntegrationMethod());

    const auto& r_N_container = this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod());

    auto nodal_temperatures = array_1d<double, TNumNodes>{};
    VariablesUtilities::GetNodalValues(this->GetGeometry(), TEMPERATURE, nodal_temperatures.begin());

    this->CalculateNodalFluxes(rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        mVariables.Np = row(r_N_container, GPoint);

        // Compute weighting coefficient for integration
        mVariables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(JContainer[GPoint],
                                                                                  IntegrationPoints[GPoint].Weight());

        this->CalculateAndAddLHS(rLeftHandSideMatrix);
        this->CalculateAndAddRHS(rRightHandSideVector, nodal_temperatures);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeProperties()
{
    KRATOS_TRY

    const auto& r_prop = this->GetProperties();

    mAlbedoCoefficient = r_prop[ALPHA_COEFFICIENT];
    mFirstCoverStorageCoefficient = r_prop[A1_COEFFICIENT];
    mSecondCoverStorageCoefficient = r_prop[A2_COEFFICIENT];
    mThirdCoverStorageCoefficient = r_prop[A3_COEFFICIENT];
    mBuildEnvironmentRadiation = r_prop[QF_COEFFICIENT];
    mMinimalStorage = r_prop[SMIN_COEFFICIENT];
    mMaximalStorage = r_prop[SMAX_COEFFICIENT];

    const GeometryType& Geom = this->GetGeometry();
    mRoughnessTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE, 1);   // This value is not read correctly, it is related to the initial value of the table
    mNetRadiation = Geom[0].FastGetSolutionStepValue(SOLAR_RADIATION, 1);  // This value is not read correctly, initial value of the table

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNetRadiation(double IncomingRadiation,
                                                                          double AtmosphericTemperature,
                                                                          double InitialSoilTemperature)
{
    // Eq 5.16
    const auto short_wave_radiation = (1.0 - mAlbedoCoefficient) * IncomingRadiation;

    constexpr auto effective_emissivity = 0.95;
    constexpr auto boltzmann_coefficient = 5.67e-8;

    // Eq 5.17
    const auto absorbed_long_wave_radiation = effective_emissivity * boltzmann_coefficient * std::pow(AtmosphericTemperature + 273.15, 4.0);

    // Eq 5.18
    const auto emitted_long_wave_radiation = boltzmann_coefficient * std::pow(InitialSoilTemperature + 273.15, 4.0);

    // Eq 5.15
    return short_wave_radiation + absorbed_long_wave_radiation - emitted_long_wave_radiation;
}

template class TMicroClimateFluxCondition<2,2>;
template class TMicroClimateFluxCondition<2,3>;
template class TMicroClimateFluxCondition<2,4>;
template class TMicroClimateFluxCondition<2,5>;
template class TMicroClimateFluxCondition<3,3>;
template class TMicroClimateFluxCondition<3,4>;
template class TMicroClimateFluxCondition<3,6>;
template class TMicroClimateFluxCondition<3,8>;
template class TMicroClimateFluxCondition<3,9>;

} // Namespace Kratos.
