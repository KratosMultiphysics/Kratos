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
#include "custom_utilities/ublas_utils.h"
#include "custom_utilities/variables_utilities.hpp"
#include "micro_climate_constants.h"

namespace Kratos
{

using namespace MicroClimateConstants;

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TMicroClimateFluxCondition<TDim, TNumNodes>::Create(
    IndexType NewId, const NodesArrayType& rNodes, Properties::Pointer pProperties) const
{
    return make_intrusive<TMicroClimateFluxCondition>(
        NewId, this->GetGeometry().Create(rNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Condition::Initialize(rCurrentProcessInfo);
    InitializeProperties();
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    CalculateRoughness(rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLocalSystem(
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
    GeometryType::JacobiansType J_container(number_of_integration_points);
    for (unsigned int i = 0; i < number_of_integration_points; ++i)
        J_container[i].resize(TDim, local_dim, false);
    r_geom.Jacobian(J_container, this->GetIntegrationMethod());

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
        const auto IntegrationCoefficient = CalculateIntegrationCoefficient(
            J_container[integration_point_index],
            r_integration_points[integration_point_index].Weight());

        CalculateAndAddLHS(rLeftHandSideMatrix, N, IntegrationCoefficient, left_hand_side_fluxes);
        CalculateAndAddRHS(rRightHandSideVector, N, IntegrationCoefficient, nodal_temperatures,
                           left_hand_side_fluxes, right_hand_side_fluxes);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
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
    mRoughnessTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE, 1);
    mNetRadiation = Geom[0].FastGetSolutionStepValue(SOLAR_RADIATION, 1);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddLHS(
    Matrix& rLeftHandSideMatrix,
    const array_1d<double, TNumNodes>& rN,
    double IntegrationCoefficient,
    const array_1d<double, TNumNodes>& rLeftHandSideFluxes)
{
    KRATOS_TRY

    auto temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{
        outer_prod(rN, rN) * IntegrationCoefficient};

    const auto flux_matrix = UBlasUtils::MakeDiagonalMatrix(
        rLeftHandSideFluxes.begin(), rLeftHandSideFluxes.end());

    temporary_matrix = prod(temporary_matrix, flux_matrix);

    GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, temporary_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(
    Vector& rRightHandSideVector,
    const array_1d<double, TNumNodes>& rN,
    double IntegrationCoefficient,
    const Vector& rNodalTemperatures,
    const array_1d<double, TNumNodes>& rLeftHandSideFluxes,
    const array_1d<double, TNumNodes>& rRightHandSideFluxes)
{
    auto temporary_matrix = BoundedMatrix<double, TNumNodes, TNumNodes>{
        outer_prod(rN, rN) * IntegrationCoefficient};
    auto temporary_vector =
        array_1d<double, TNumNodes>{prod(temporary_matrix, rRightHandSideFluxes)};
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, temporary_vector);

    const auto flux_matrix = UBlasUtils::MakeDiagonalMatrix(
        rLeftHandSideFluxes.begin(), rLeftHandSideFluxes.end());
    temporary_matrix = prod(temporary_matrix, flux_matrix);
    temporary_vector = -prod(temporary_matrix, rNodalTemperatures);
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, temporary_vector);
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const Matrix& rJacobian, double Weight) const
{
    if (TDim == 2)
    {
        return MathUtils<>::Norm(
                   UBlasUtils::MakeVector({rJacobian(0, 0), rJacobian(1, 0)})) *
               Weight;
    }
    else if (TDim == 3)
    {
        const auto vec1 =
            UBlasUtils::MakeVector({rJacobian(0, 0), rJacobian(1, 0), rJacobian(2, 0)});
        const auto vec2 =
            UBlasUtils::MakeVector({rJacobian(0, 1), rJacobian(1, 1), rJacobian(2, 1)});
        return MathUtils<>::Norm(MathUtils<>::CrossProduct(vec1, vec2)) * Weight;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLeftHandSideFluxes() const
{
    array_1d<double, TNumNodes> result;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        // Eq 5.22
        const auto sensible_heat_flux_left =
            AirHeatCapacity * AirDensity / RoughnessLayerResistance;
        result[i] = sensible_heat_flux_left;
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRightHandSideFluxes(
    double time_step_size, double previous_storage, double previous_radiation) const
{
    array_1d<double, TNumNodes> result;
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto net_radiation = CalculateNetRadiation(i);
        const auto surface_heat_storage = CalculateSurfaceHeatStorage(
            time_step_size, previous_radiation, net_radiation);
        const WaterFluxes water_fluxes = CalculateWaterFluxes(
            i, time_step_size, previous_storage, net_radiation, surface_heat_storage);

        result[i] = CalculateRightHandSideFlux(
            net_radiation, surface_heat_storage, water_fluxes.evaporation);
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRightHandSideFlux(
    double net_radiation, double surface_heat_storage, double actual_evaporation) const
{
    const auto latent_heat_flux = actual_evaporation * WaterDensity * LatentEvaporationHeat;

    // Eq 5.22
    const auto sensible_heat_flux_right =
        -AirHeatCapacity * AirDensity * mRoughnessTemperature / RoughnessLayerResistance;

    // Eq 5.31
    const auto subsurface_heat_flux =
        net_radiation - sensible_heat_flux_right - latent_heat_flux +
        mBuildEnvironmentRadiation - surface_heat_storage;
    return subsurface_heat_flux;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateCurrentWaterStorage(
    double time_step_size, double previous_storage, double previous_radiation) const
{
    double result = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto net_radiation = CalculateNetRadiation(i);
        const auto surface_heat_storage = CalculateSurfaceHeatStorage(
            time_step_size, previous_radiation, net_radiation);
        const WaterFluxes water_fluxes = CalculateWaterFluxes(
            i, time_step_size, previous_storage, net_radiation, surface_heat_storage);

        const auto actual_storage =
            previous_storage +
            time_step_size * (water_fluxes.precipitation - water_fluxes.evaporation);
        result += actual_storage / TNumNodes;
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateCurrentNetRadiation() const
{
    double result = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto net_radiation = CalculateNetRadiation(i);
        result += net_radiation / TNumNodes;
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNetRadiation(unsigned int index) const
{
    auto& r_geom = this->GetGeometry();

    // Eq 5.16
    const auto incoming_radiation = r_geom[index].FastGetSolutionStepValue(SOLAR_RADIATION);
    const auto short_wave_radiation = (1.0 - mAlbedoCoefficient) * incoming_radiation;

    // Eq 5.17
    const auto atmospheric_temperature =
        r_geom[index].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto absorbed_long_wave_radiation =
        EffectiveEmissivity * BoltzmannCoefficient *
        std::pow(atmospheric_temperature + 273.15, 4.0);

    // Eq 5.18
    const auto initial_soil_temperature =
        r_geom[index].FastGetSolutionStepValue(TEMPERATURE, 1);
    const auto emitted_long_wave_radiation =
        BoltzmannCoefficient * std::pow(initial_soil_temperature + 273.15, 4.0);

    // Eq 5.15
    return short_wave_radiation + absorbed_long_wave_radiation - emitted_long_wave_radiation;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateSurfaceHeatStorage(
    double time_step_size, double previous_radiation, double net_radiation) const
{
    // Eq 5.20
    return mFirstCoverStorageCoefficient * net_radiation +
           mSecondCoverStorageCoefficient * (net_radiation - previous_radiation) / time_step_size +
           mThirdCoverStorageCoefficient;
}

template <unsigned int TDim, unsigned int TNumNodes>
WaterFluxes TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateWaterFluxes(
    unsigned int i, double time_step_size, double previous_storage, double net_radiation, double surface_heat_storage) const
{
    const auto potential_evaporation =
        CalculatePotentialEvaporation(i, net_radiation, surface_heat_storage);

    auto& r_geom = this->GetGeometry();
    const auto precipitation = r_geom[i].FastGetSolutionStepValue(PRECIPITATION);
    const auto potential_storage =
        previous_storage + time_step_size * (precipitation - potential_evaporation);

    auto actual_precipitation = precipitation;
    auto actual_evaporation = potential_evaporation;

    // Correct the precipitation and evaporation if the storage is out of bounds
    if (potential_storage > mMaximalStorage)
    {
        actual_precipitation =
            (mMaximalStorage - previous_storage) / time_step_size + actual_evaporation;
    }
    else if (potential_storage < mMinimalStorage)
    {
        actual_evaporation =
            (previous_storage - mMinimalStorage) / time_step_size + actual_precipitation;
    }

    WaterFluxes water_fluxes;
    water_fluxes.evaporation = actual_evaporation;
    water_fluxes.precipitation = actual_precipitation;
    return water_fluxes;
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculatePotentialEvaporation(
    unsigned int i, double net_radiation, double surface_heat_storage) const
{
    auto& r_geom = this->GetGeometry();

    // Eq 5.35
    const auto wind_speed = r_geom[i].FastGetSolutionStepValue(WIND_SPEED);
    const auto atmospheric_resistance = 1.0 / (0.007 + 0.0056 * wind_speed);

    // Eq. 5.12
    const auto atmospheric_temperature = r_geom[i].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto saturated_vapor_pressure =
        6.11 * std::exp(17.27 * atmospheric_temperature / (atmospheric_temperature + 237.3));

    // Eq 5.13
    const auto vapor_pressure_increment =
        4098.0 * saturated_vapor_pressure /
        (std::pow((atmospheric_temperature + 237.3), 2.0));

    // Eq 5.14
    const auto humidity = r_geom[i].FastGetSolutionStepValue(AIR_HUMIDITY);
    const auto actual_vapor_pressure = humidity / 100.0 * saturated_vapor_pressure;

    // Eq 5.34
    auto latent_heat_flux =
        (vapor_pressure_increment * (net_radiation + mBuildEnvironmentRadiation - surface_heat_storage) +
         AirHeatCapacity * AirDensity *
             (saturated_vapor_pressure - actual_vapor_pressure) / atmospheric_resistance) /
        (vapor_pressure_increment +
         PsychometricConstant * (1.0 + SurfaceResistance / atmospheric_resistance));
    latent_heat_flux = std::max(latent_heat_flux, 0.0);

    // Eq 5.36
    return latent_heat_flux / (WaterDensity * LatentEvaporationHeat);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRoughness(const ProcessInfo& rCurrentProcessInfo)
{
    const auto time_step_size = rCurrentProcessInfo.GetValue(DELTA_TIME);

    const auto& r_geom = this->GetGeometry();
    const auto current_air_temperature = r_geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE);
    const auto current_wind_speed =
        std::max(1.0e-3, r_geom[0].FastGetSolutionStepValue(WIND_SPEED));

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

        const auto c = RoughnessLayerResistance * RoughnessLayerHeight + time_step_size +
                       time_step_size * current_wind_speed *
                           RoughnessLayerResistance * surface_roughness_factor *
                           friction_drag_coefficient * friction_drag_coefficient;
        const auto current_roughness_temperature =
            (RoughnessLayerResistance * RoughnessLayerHeight * previous_roughness_temperature +
             time_step_size * initial_soil_temperature +
             time_step_size * current_wind_speed * RoughnessLayerResistance *
                 surface_roughness_factor * friction_drag_coefficient *
                 friction_drag_coefficient * current_air_temperature) /
            c;
        mRoughnessTemperature += current_roughness_temperature / TNumNodes;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateSurfaceRoughnessFactor(
    double CurrentAirTemperature,
    double PreviousRoughnessTemperature,
    double RichardsonBulkModulus,
    double FrictionDragCoefficient) const
{
    if (PreviousRoughnessTemperature >= CurrentAirTemperature)
    {
        // Eq 5.27
        const auto cof = RichardsonBulkModulus /
                         (1.0 + 75.0 * FrictionDragCoefficient * FrictionDragCoefficient *
                                    std::sqrt(MeasurementHeight / RoughnessHeight *
                                              std::abs(RichardsonBulkModulus)));
        return 1.0 - 15.0 * cof;
    }

    // Eq 5.28
    const auto cof = std::sqrt(1.0 + 5.0 * RichardsonBulkModulus);
    return 1.0 / (1.0 + 15.0 * RichardsonBulkModulus * cof);
}

template class TMicroClimateFluxCondition<2, 2>;
template class TMicroClimateFluxCondition<2, 3>;
template class TMicroClimateFluxCondition<2, 4>;
template class TMicroClimateFluxCondition<2, 5>;
template class TMicroClimateFluxCondition<3, 3>;
template class TMicroClimateFluxCondition<3, 4>;
template class TMicroClimateFluxCondition<3, 6>;
template class TMicroClimateFluxCondition<3, 8>;
template class TMicroClimateFluxCondition<3, 9>;

} // Namespace Kratos.
