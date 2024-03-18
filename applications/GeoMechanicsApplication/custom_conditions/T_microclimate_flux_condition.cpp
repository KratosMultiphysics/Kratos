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
#include "custom_conditions/T_microclimate_flux_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "micro_climate_constants.h"

#include <boost/numeric/ublas/vector_expression.hpp>

namespace Kratos
{

using namespace MicroClimateConstants;

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer GeoTMicroClimateFluxCondition<TDim, TNumNodes>::Create(
    IndexType NewId, const NodesArrayType& rNodes, Properties::Pointer pProperties) const
{
    return make_intrusive<GeoTMicroClimateFluxCondition>(
        NewId, this->GetGeometry().Create(rNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    GeoTCondition<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);
    InitializeProperties();
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // Unfortunately, the roughness temperature and net radiation cannot be
    // initialized in the Initialize() method because the values of the air
    // temperature and solar radiation are not yet available.
    if (!mIsInitialized)
    {
        const GeometryType& Geom = this->GetGeometry();
        mRoughnessTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE, 0);
        mNetRadiation = Geom[0].FastGetSolutionStepValue(SOLAR_RADIATION, 0);
        mIsInitialized = true;
    }

    CalculateRoughness(rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLocalSystem(
    Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rLeftHandSideMatrix = Matrix{TNumNodes, TNumNodes, 0.0};
    rRightHandSideVector = Vector{TNumNodes, 0.0};

    // Previous definitions
    const auto& r_geom = this->GetGeometry();
    const auto& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const auto number_of_integration_points =
        static_cast<unsigned int>(r_integration_points.size());

    // Containers of variables at all integration points
    const auto local_dim = static_cast<unsigned int>(r_geom.LocalSpaceDimension());
    GeometryType::JacobiansType jacobians(number_of_integration_points);
    for (unsigned int i = 0; i < number_of_integration_points; ++i)
        jacobians[i].resize(TDim, local_dim, false);
    r_geom.Jacobian(jacobians, this->GetIntegrationMethod());

    const auto& r_N_container =
        this->GetGeometry().ShapeFunctionsValues(this->GetIntegrationMethod());

    auto nodal_temperatures = array_1d<double, TNumNodes>{};
    VariablesUtilities::GetNodalValues(this->GetGeometry(), TEMPERATURE,
                                       nodal_temperatures.begin());

    const auto time_step_size = rCurrentProcessInfo.GetValue(DELTA_TIME);

    const auto previous_storage = mWaterStorage;
    const auto previous_radiation = mNetRadiation;
    mWaterStorage = CalculateCurrentWaterStorage(
        time_step_size, previous_storage, previous_radiation);
    mNetRadiation = CalculateCurrentNetRadiation();

    const auto left_hand_side_fluxes = CalculateLeftHandSideFluxes();
    const auto right_hand_side_fluxes = CalculateRightHandSideFluxes(
        time_step_size, previous_storage, previous_radiation);

    // Loop over integration points
    for (unsigned int integration_point_index = 0;
         integration_point_index < number_of_integration_points; ++integration_point_index)
    {
        const auto N =
            array_1d<double, TNumNodes>{row(r_N_container, integration_point_index)};

        // Compute weighting coefficient for integration
        const auto IntegrationCoefficient =
            ConditionUtilities::CalculateIntegrationCoefficient<TDim, TNumNodes>(
                jacobians[integration_point_index],
                r_integration_points[integration_point_index].Weight());

        CalculateAndAddLHS(rLeftHandSideMatrix, N, IntegrationCoefficient, left_hand_side_fluxes);
        CalculateAndAddRHS(rRightHandSideVector, N, IntegrationCoefficient, nodal_temperatures,
                           left_hand_side_fluxes, right_hand_side_fluxes);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::InitializeProperties()
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
    mWaterDensity = r_prop[DENSITY_WATER];

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddLHS(
    Matrix& rLeftHandSideMatrix,
    const array_1d<double, TNumNodes>& rN,
    double IntegrationCoefficient,
    const array_1d<double, TNumNodes>& rLeftHandSideFluxes)
{
    KRATOS_TRY

    const auto flux_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{
        outer_prod(rN, element_prod(rN, rLeftHandSideFluxes)) * IntegrationCoefficient};

    rLeftHandSideMatrix += flux_matrix;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(
    Vector& rRightHandSideVector,
    const array_1d<double, TNumNodes>& rN,
    double IntegrationCoefficient,
    const Vector& rNodalTemperatures,
    const array_1d<double, TNumNodes>& rLeftHandSideFluxes,
    const array_1d<double, TNumNodes>& rRightHandSideFluxes)
{
    auto temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{
        outer_prod(rN, rN) * IntegrationCoefficient};
    rRightHandSideVector += prod(temporary_matrix, rRightHandSideFluxes);

    temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{
        outer_prod(rN, element_prod(rN, rLeftHandSideFluxes)) * IntegrationCoefficient};
    rRightHandSideVector -= prod(temporary_matrix, rNodalTemperatures);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLeftHandSideFluxes() const
{
    // Eq 5.22
    const auto sensible_heat_flux_left = AirHeatCapacity * AirDensity / RoughnessLayerResistance;
    return array_1d<double, TNumNodes>(TNumNodes, sensible_heat_flux_left);
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRightHandSideFluxes(
    double TimeStepSize, double PreviousStorage, double PreviousRadiation) const
{
    array_1d<double, TNumNodes> result;
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto net_radiation = CalculateNetRadiation(i);
        const auto surface_heat_storage = CalculateSurfaceHeatStorage(
            TimeStepSize, PreviousRadiation, net_radiation);
        const WaterFluxes water_fluxes = CalculateWaterFluxes(
            i, TimeStepSize, PreviousStorage, net_radiation, surface_heat_storage);

        result[i] = CalculateRightHandSideFlux(
            net_radiation, surface_heat_storage, water_fluxes.evaporation);
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRightHandSideFlux(
    double NetRadiation, double SurfaceHeatStorage, double ActualEvaporation) const
{
    const auto latent_heat_flux = ActualEvaporation * mWaterDensity * LatentEvaporationHeat;

    // Eq 5.22
    const auto sensible_heat_flux_right =
        -AirHeatCapacity * AirDensity * mRoughnessTemperature / RoughnessLayerResistance;

    // Eq 5.31
    const auto subsurface_heat_flux =
        NetRadiation - sensible_heat_flux_right - latent_heat_flux +
        mBuildEnvironmentRadiation - SurfaceHeatStorage;
    return subsurface_heat_flux;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateCurrentWaterStorage(
    double TimeStepSize, double PreviousStorage, double PreviousRadiation) const
{
    double result = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto net_radiation = CalculateNetRadiation(i);
        const auto surface_heat_storage = CalculateSurfaceHeatStorage(
            TimeStepSize, PreviousRadiation, net_radiation);
        const WaterFluxes water_fluxes = CalculateWaterFluxes(
            i, TimeStepSize, PreviousStorage, net_radiation, surface_heat_storage);

        const auto actual_storage =
            PreviousStorage +
            TimeStepSize * (water_fluxes.precipitation - water_fluxes.evaporation);
        result += actual_storage;
    }

    return result / TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateCurrentNetRadiation() const
{
    std::vector<unsigned int> local_node_indices(TNumNodes);
    std::iota(local_node_indices.begin(), local_node_indices.end(), 0u);
    const auto result =
        std::accumulate(local_node_indices.begin(), local_node_indices.end(), 0.0,
                        [this](double sum, unsigned int i)
                        { return sum + CalculateNetRadiation(i); }) /
        TNumNodes;
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNetRadiation(unsigned int NodeIndex) const
{
    auto& r_geom = this->GetGeometry();

    // Eq 5.16
    const auto incoming_radiation =
        r_geom[NodeIndex].FastGetSolutionStepValue(SOLAR_RADIATION);
    const auto short_wave_radiation = (1.0 - mAlbedoCoefficient) * incoming_radiation;

    // Eq 5.17
    const auto atmospheric_temperature =
        r_geom[NodeIndex].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto absorbed_long_wave_radiation =
        EffectiveEmissivity * BoltzmannCoefficient *
        std::pow(atmospheric_temperature + 273.15, 4.0);

    // Eq 5.18
    const auto initial_soil_temperature =
        r_geom[NodeIndex].FastGetSolutionStepValue(TEMPERATURE, 1);
    const auto emitted_long_wave_radiation =
        BoltzmannCoefficient * std::pow(initial_soil_temperature + 273.15, 4.0);

    // Eq 5.15
    return short_wave_radiation + absorbed_long_wave_radiation - emitted_long_wave_radiation;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateSurfaceHeatStorage(
    double TimeStepSize, double PreviousRadiation, double NetRadiation) const
{
    // Eq 5.20
    return mFirstCoverStorageCoefficient * NetRadiation +
           mSecondCoverStorageCoefficient * (NetRadiation - PreviousRadiation) / TimeStepSize +
           mThirdCoverStorageCoefficient;
}

template <unsigned int TDim, unsigned int TNumNodes>
WaterFluxes GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateWaterFluxes(
    unsigned int NodeIndex, double TimeStepSize, double PreviousStorage, double NetRadiation, double SurfaceHeatStorage) const
{
    const auto potential_evaporation =
        CalculatePotentialEvaporation(NodeIndex, NetRadiation, SurfaceHeatStorage);

    auto& r_geom = this->GetGeometry();
    const auto precipitation = r_geom[NodeIndex].FastGetSolutionStepValue(PRECIPITATION);
    const auto potential_storage =
        PreviousStorage + TimeStepSize * (precipitation - potential_evaporation);

    auto actual_precipitation = precipitation;
    auto actual_evaporation = potential_evaporation;

    // Correct the precipitation and evaporation if the storage is out of bounds
    if (potential_storage > mMaximalStorage)
    {
        actual_precipitation =
            (mMaximalStorage - PreviousStorage) / TimeStepSize + actual_evaporation;
    }
    else if (potential_storage < mMinimalStorage)
    {
        actual_evaporation = (PreviousStorage - mMinimalStorage) / TimeStepSize + actual_precipitation;
    }

    return {actual_precipitation, actual_evaporation};
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculatePotentialEvaporation(
    unsigned int NodeIndex, double NetRadiation, double SurfaceHeatStorage) const
{
    auto& r_geom = this->GetGeometry();

    // Eq 5.35
    const auto wind_speed = r_geom[NodeIndex].FastGetSolutionStepValue(WIND_SPEED);
    const auto atmospheric_resistance = 1.0 / (0.007 + 0.0056 * wind_speed);

    // Eq. 5.12
    const auto atmospheric_temperature =
        r_geom[NodeIndex].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto saturated_vapor_pressure =
        6.11 * std::exp(17.27 * atmospheric_temperature / (atmospheric_temperature + 237.3));

    // Eq 5.13
    const auto vapor_pressure_increment =
        4098.0 * saturated_vapor_pressure /
        (std::pow((atmospheric_temperature + 237.3), 2.0));

    // Eq 5.14
    const auto humidity = r_geom[NodeIndex].FastGetSolutionStepValue(AIR_HUMIDITY);
    const auto actual_vapor_pressure = (humidity / 100.0) * saturated_vapor_pressure;

    // Eq 5.34
    auto latent_heat_flux =
        (vapor_pressure_increment * (NetRadiation + mBuildEnvironmentRadiation - SurfaceHeatStorage) +
         AirHeatCapacity * AirDensity *
             (saturated_vapor_pressure - actual_vapor_pressure) / atmospheric_resistance) /
        (vapor_pressure_increment +
         PsychometricConstant * (1.0 + SurfaceResistance / atmospheric_resistance));
    latent_heat_flux = std::max(latent_heat_flux, 0.0);

    // Eq 5.36
    return latent_heat_flux / (mWaterDensity * LatentEvaporationHeat);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRoughness(const ProcessInfo& rCurrentProcessInfo)
{
    const auto time_step_size = rCurrentProcessInfo.GetValue(DELTA_TIME);

    const auto& r_geom = this->GetGeometry();
    const auto current_air_temperature = r_geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto lower_bound = 1.0e-3;
    const auto current_wind_speed =
        std::max(lower_bound, r_geom[0].FastGetSolutionStepValue(WIND_SPEED));

    const auto previous_roughness_temperature = mRoughnessTemperature;
    mRoughnessTemperature = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto initial_soil_temperature =
            r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

        // Eq 5.29
        const auto richardson_bulk_modulus =
            2.0 * GravitationalAcceleration * MeasurementHeight /
            (current_air_temperature + previous_roughness_temperature + 546.3) *
            (current_air_temperature - previous_roughness_temperature) /
            (current_wind_speed * current_wind_speed);

        // Eq 5.25
        const auto friction_drag_coefficient =
            VonNeumannCoefficient / std::log(MeasurementHeight / RoughnessHeight);

        const auto surface_roughness_factor = CalculateSurfaceRoughnessFactor(
            current_air_temperature, previous_roughness_temperature,
            richardson_bulk_modulus, friction_drag_coefficient);

        // Eq 5.30
        const auto denominator =
            RoughnessLayerResistance * RoughnessLayerHeight + time_step_size +
            time_step_size * current_wind_speed * RoughnessLayerResistance *
                surface_roughness_factor * friction_drag_coefficient * friction_drag_coefficient;
        const auto current_roughness_temperature =
            (RoughnessLayerResistance * RoughnessLayerHeight * previous_roughness_temperature +
             time_step_size * initial_soil_temperature +
             time_step_size * current_wind_speed * RoughnessLayerResistance *
                 surface_roughness_factor * friction_drag_coefficient *
                 friction_drag_coefficient * current_air_temperature) /
            denominator;
        mRoughnessTemperature += current_roughness_temperature;
    }

    mRoughnessTemperature /= TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoTMicroClimateFluxCondition<TDim, TNumNodes>::CalculateSurfaceRoughnessFactor(
    double CurrentAirTemperature,
    double PreviousRoughnessTemperature,
    double RichardsonBulkModulus,
    double FrictionDragCoefficient) const
{
    if (PreviousRoughnessTemperature >= CurrentAirTemperature)
    {
        // Eq 5.27
        const auto coefficient =
            RichardsonBulkModulus /
            (1.0 + 75.0 * FrictionDragCoefficient * FrictionDragCoefficient *
                       std::sqrt(MeasurementHeight / RoughnessHeight *
                                 std::abs(RichardsonBulkModulus)));
        return 1.0 - 15.0 * coefficient;
    }

    // Eq 5.28
    const auto coefficient = std::sqrt(1.0 + 5.0 * RichardsonBulkModulus);
    return 1.0 / (1.0 + 15.0 * RichardsonBulkModulus * coefficient);
}

template class GeoTMicroClimateFluxCondition<2, 2>;
template class GeoTMicroClimateFluxCondition<2, 3>;
template class GeoTMicroClimateFluxCondition<2, 4>;
template class GeoTMicroClimateFluxCondition<2, 5>;
template class GeoTMicroClimateFluxCondition<3, 3>;
template class GeoTMicroClimateFluxCondition<3, 4>;
template class GeoTMicroClimateFluxCondition<3, 6>;
template class GeoTMicroClimateFluxCondition<3, 8>;
template class GeoTMicroClimateFluxCondition<3, 9>;

} // Namespace Kratos.
