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
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAll(
    Matrix& rLeftHandSideMatrix,
    Vector& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
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
    const ProcessInfo& CurrentProcessInfo)
{
    double timeStepSize = CurrentProcessInfo.GetValue(DELTA_TIME);

    const GeometryType& Geom = this->GetGeometry();
    const double currentAirTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE);
    double currentWindSpeed = Geom[0].FastGetSolutionStepValue(WIND_SPEED);

    constexpr double roughnessLayerHeight = 10.0;
    constexpr double roughnessLayerResistance = 30.0;
    constexpr double vonNeumanCoefficient = 0.4;
    constexpr double measurementHeight = 10.0;
    constexpr double roughnessHeight = 1.0;
    constexpr double gravitationalAcceleration = 9.81;

    const auto previous_roughness_temperature = mVariables.roughnessTemperature;
    mVariables.roughnessTemperature = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const double initialSoilTemperature = Geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

        double surfaceRoughnessFactor = 0.0;

        currentWindSpeed = std::max(currentWindSpeed, 1.0e-3);

        // Eq 5.29
        const double richardsonBulkModulus = 2.0 * gravitationalAcceleration * measurementHeight / (currentAirTemperature +
            previous_roughness_temperature + 546.3) * (currentAirTemperature - previous_roughness_temperature) / (currentWindSpeed * currentWindSpeed);

        // Eq 5.25
        const double frictionDragCoefficient = vonNeumanCoefficient / std::log(measurementHeight / roughnessHeight);

        double cof = 0.0;
        if (previous_roughness_temperature >= currentAirTemperature) {
            // Eq 5.27
            cof = richardsonBulkModulus / (1.0 + 75.0 * frictionDragCoefficient * frictionDragCoefficient *
                std::sqrt(measurementHeight / roughnessHeight * std::fabs(richardsonBulkModulus)));
            surfaceRoughnessFactor = 1.0 - 15.0 * cof;
        }
        else {
            // Eq 5.28
            cof = std::sqrt(1.0 + 5.0 * richardsonBulkModulus);
            surfaceRoughnessFactor = 1.0 / (1.0 + 15.0 * richardsonBulkModulus * cof);
        }

        const double c = roughnessLayerResistance * roughnessLayerHeight + timeStepSize + timeStepSize * currentWindSpeed * roughnessLayerResistance *
            surfaceRoughnessFactor * frictionDragCoefficient * frictionDragCoefficient;
        double currentRoughnessTemperature = (roughnessLayerResistance * roughnessLayerHeight * previous_roughness_temperature + timeStepSize *
            initialSoilTemperature + timeStepSize * currentWindSpeed * roughnessLayerResistance * surfaceRoughnessFactor *
            frictionDragCoefficient * frictionDragCoefficient * currentAirTemperature) / c;
        mVariables.roughnessTemperature += currentRoughnessTemperature / TNumNodes;
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNodalFluxes(
    const ProcessInfo& CurrentProcessInfo)
{
    const double maximalStorage = mVariables.maximalStorage;

    const Properties mProperties = this->GetProperties();

    const GeometryType& Geom = this->GetGeometry();
    const double timeStepSize = CurrentProcessInfo.GetValue(DELTA_TIME);

    constexpr double airDensity = 1.18;
    constexpr double airHeatCapacity = 1004.67;
    constexpr double roughnessLayerResistance = 30.0;
    constexpr double latentEvaporationHeat = 2.45e6;
    constexpr double boltzmannCoefficient = 5.67e-8;
    constexpr double effectiveEmissivity = 0.95;
    constexpr double waterDensity = 1e3;
    constexpr double psychometricConstant = 0.63;
    constexpr double surfaceResistance = 30.0;

    const auto previous_storage = mVariables.waterStorage;
    const auto previous_radiation = mVariables.netRadiation;
    mVariables.waterStorage = 0.0;
    mVariables.netRadiation = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const double atmosphericTemperature = Geom[i].FastGetSolutionStepValue(AIR_TEMPERATURE);
        const double incomingRadiation = Geom[i].FastGetSolutionStepValue(SOLAR_RADIATION);
        const double humidity = Geom[i].FastGetSolutionStepValue(AIR_HUMIDITY);
        const double precipitation = Geom[i].FastGetSolutionStepValue(PRECIPITATION);
        const double windSpeed = Geom[i].FastGetSolutionStepValue(WIND_SPEED);

        const double initialSoilTemperature = Geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

        // Eq 5.22
        const double sensibleHeatFluxLeft = airHeatCapacity * airDensity / roughnessLayerResistance;
        const double sensibleHeatFluxRight = -airHeatCapacity * airDensity * mVariables.roughnessTemperature / roughnessLayerResistance;

        // Eq 5.35
        const double atmosphericResistance = 1.0 / (0.007 + 0.0056 * windSpeed);

        // Eq. 5.12
        const double saturatedVaporPressure = 6.11 * std::exp(17.27 * atmosphericTemperature / (atmosphericTemperature + 237.3));

        // Eq 5.13
        const double vaporPressureIncrement = 4098.0 * saturatedVaporPressure / (std::pow((atmosphericTemperature + 237.3), 2.0));

        // Eq 5.14
        const double actualVaporPressure = humidity / 100.0 * saturatedVaporPressure;

        // Eq 5.16
        const double shortWaveRadiation = (1.0 - mAlbedoCoefficient) * incomingRadiation;

        // Eq 5.18
        const double emittedLongWaveRadiation = boltzmannCoefficient * std::pow(initialSoilTemperature + 273.15, 4.0);

        // Eq 5.17
        const double absorbedLongWaveRadiation = effectiveEmissivity * boltzmannCoefficient * std::pow(atmosphericTemperature + 273.15, 4.0);

        // Eq 5.15
        const double netRadiation = shortWaveRadiation + absorbedLongWaveRadiation - emittedLongWaveRadiation;

        // Eq 5.20
        const double surfaceHeatStorage = mFirstCoverStorageCoefficient * netRadiation + mSecondCoverStorageCoefficient * (netRadiation - previous_radiation) /
            timeStepSize + mThirdCoverStorageCoefficient;

        // Eq 5.34
        double latentHeatFlux = (vaporPressureIncrement * (netRadiation + mBuildEnvironmentRadiation - surfaceHeatStorage) + airHeatCapacity * airDensity *
            (saturatedVaporPressure - actualVaporPressure) / atmosphericResistance) / (vaporPressureIncrement + psychometricConstant *  // division is replaced to * based on (3.34)
                (1.0 + surfaceResistance / atmosphericResistance));   //Where is G (5.34)?
        latentHeatFlux = std::max(latentHeatFlux, 0.0);

        double actualEvaporation = 0.0;
        double actualPrecipitation = 0.0;

        // Eq 5.36
        const double potentialEvaporation = latentHeatFlux / (waterDensity * latentEvaporationHeat);
        const double potentialStorage = previous_storage + timeStepSize * (precipitation - potentialEvaporation);
        if (potentialStorage > maximalStorage)
        {
            actualEvaporation = potentialEvaporation;
            actualPrecipitation = (maximalStorage - previous_storage) / timeStepSize + actualEvaporation;
        }
        else if (potentialStorage < mMinimalStorage)
        {
            actualPrecipitation = precipitation;
            actualEvaporation = (previous_storage - mMinimalStorage) / timeStepSize + actualPrecipitation;
        }
        else
        {
            actualEvaporation = potentialEvaporation;
            actualPrecipitation = precipitation;
        }
        const double actualStorage = previous_storage + timeStepSize * (actualPrecipitation - actualEvaporation);
        latentHeatFlux = actualEvaporation * waterDensity * latentEvaporationHeat;

        // Eq 5.31
        double subsurfaceHeatFlux = netRadiation - sensibleHeatFluxRight - latentHeatFlux + mBuildEnvironmentRadiation - surfaceHeatStorage;

        mVariables.netRadiation += netRadiation / TNumNodes;
        mVariables.waterStorage += actualStorage / TNumNodes;
        mVariables.leftHandSideFlux[i] = sensibleHeatFluxLeft;
        mVariables.rightHandSideFlux[i] = subsurfaceHeatFlux;
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
    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

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
    mVariables.maximalStorage = r_prop[SMAX_COEFFICIENT];

    const GeometryType& Geom = this->GetGeometry();
    mVariables.roughnessTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE, 1);   // This value is not read correctly
    mVariables.waterStorage = 0.0;                                                            // it is related to the initial value of the table
    mVariables.netRadiation = Geom[0].FastGetSolutionStepValue(SOLAR_RADIATION, 1);  // This value is not read correctly, initial value of the table

    KRATOS_CATCH("")
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
