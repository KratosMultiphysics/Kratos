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
//


// Application includes
#include "custom_conditions/T_microclimate_flux_condition.hpp"

namespace Kratos
{

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TMicroClimateFluxCondition<TDim,TNumNodes>::Create(
    IndexType NewId,NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new TMicroClimateFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    if (!mIsInitialised)
    {
        this->InitializeProperties();
        mIsInitialised = true;
    }
    this->CalculateRoughness(rCurrentProcessInfo, rVariables);
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Containers of variables at all integration points
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    Geom.Jacobian(JContainer, this->GetIntegrationMethod());

    //Element variables
    this->InitializeElementVariables(rVariables, rCurrentProcessInfo);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(rVariables, GPoint);

        //Compute weighting coefficient for integration
        //Variables.IntegrationCoefficient =
        //    this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);
        this->CalculateIntegrationCoefficient(rVariables.IntegrationCoefficient,
            JContainer[GPoint],
            IntegrationPoints[GPoint].Weight());

        //Contributions to the left hand side
        this->CalculateAndAddLHS(rLeftHandSideMatrix, rVariables);

        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, rVariables);
    }

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeElementVariables(
    ElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Properties variables
    //this->InitializeProperties(rVariables);

    //Nodal Variables
    this->InitializeNodalTemperatureVariables(rVariables);

    //Nodal Variables
    //this->CalculateRoughness(rCurrentProcessInfo, rVariables);
    this->CalculateNodalFluxes(rCurrentProcessInfo, rVariables);

    //Variables computed at each GP
    rVariables.Np.resize(TNumNodes, false);

    //General Variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    // shape functions
    (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
    rVariables.NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    // gradient of shape functions and determinant of Jacobian
    rVariables.detJContainer.resize(NumGPoints, false);

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeNodalTemperatureVariables(
    ElementVariables& rVariables)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    //Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rVariables.TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        rVariables.DtTemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_TEMPERATURE);
    }

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    rVariables.TMatrix = outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;
    rVariables.TVector = prod(rVariables.TMatrix, rVariables.rightHandSideFlux);
    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.TVector);

    //---------------

    rVariables.TMatrix = outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;
    Matrix TTMatrix = ZeroMatrix(TNumNodes, TNumNodes);
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        TTMatrix(i, i) = rVariables.leftHandSideFlux[i];
    }
    rVariables.TMatrix = prod(rVariables.TMatrix, TTMatrix);
    rVariables.TVector = -prod(rVariables.TMatrix, rVariables.TemperatureVector);
    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.TVector);
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    ElementVariables& rVariables)
{
    KRATOS_TRY

    rVariables.TMatrix = outer_prod(rVariables.Np, rVariables.Np) * rVariables.IntegrationCoefficient;

    Matrix TTMatrix = ZeroMatrix(TNumNodes, TNumNodes);
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        TTMatrix(i, i) = rVariables.leftHandSideFlux[i];
    }

    rVariables.TMatrix = prod(rVariables.TMatrix, TTMatrix);

    GeoElementUtilities::
        AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rVariables.TMatrix);

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim,TNumNodes>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    if (TDim == 2)
    {
        const double dx_dxi = Jacobian(0, 0);
        const double dy_dxi = Jacobian(1, 0);
        const double ds = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);
        rIntegrationCoefficient = ds * Weight;
    }
    else if (TDim == 3)
    {
        double NormalVector[3];
        NormalVector[0] = Jacobian(1, 0) * Jacobian(2, 1) - Jacobian(2, 0) * Jacobian(1, 1);
        NormalVector[1] = Jacobian(2, 0) * Jacobian(0, 1) - Jacobian(0, 0) * Jacobian(2, 1);
        NormalVector[2] = Jacobian(0, 0) * Jacobian(1, 1) - Jacobian(1, 0) * Jacobian(0, 1);
        const double dA = std::sqrt(NormalVector[0] * NormalVector[0]
            + NormalVector[1] * NormalVector[1]
            + NormalVector[2] * NormalVector[2]);
        rIntegrationCoefficient = dA * Weight;
    }
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateRoughness(
    const ProcessInfo& CurrentProcessInfo,
    ElementVariables& rVariables)
{
    double timeStepSize = CurrentProcessInfo.GetValue(DELTA_TIME);

    const GeometryType& Geom = this->GetGeometry();
    const double currentAirTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE);
    double currentWindSpeed = Geom[0].FastGetSolutionStepValue(WIND_SPEED);

    const Properties mProperties = this->GetProperties();

    constexpr double roughnessLayerHeight = 10.0;
    constexpr double roughnessLayerResistance = 30.0;
    constexpr double vonNeumanCoefficient = 0.4;
    constexpr double measurementHeight = 10.0;
    constexpr double roughnessHeight = 1.0;
    constexpr double gravitationalAcceleration = 9.81;

    rVariables.previousRoughnessTemperature = rVariables.roughnessTemperature;
    rVariables.previousStorage = rVariables.waterStorage;
    rVariables.previousRadiation = rVariables.netRadiation;
    rVariables.roughnessTemperature = 0.0;

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const double initialSoilTemperature = Geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);

        double surfaceRoughnessFactor = 0.0;

        currentWindSpeed = std::max(currentWindSpeed, 1.0e-3);

        // Eq 5.29
        const double richardsonBulkModulus = 2.0 * gravitationalAcceleration * measurementHeight / (currentAirTemperature +
            rVariables.previousRoughnessTemperature + 546.3) * (currentAirTemperature - rVariables.previousRoughnessTemperature) / (currentWindSpeed * currentWindSpeed);

        // Eq 5.25
        const double frictionDragCoefficient = vonNeumanCoefficient / std::log(measurementHeight / roughnessHeight);

        double cof = 0.0;
        if (rVariables.previousRoughnessTemperature >= currentAirTemperature) {
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
        double currentRoughnessTemperature = (roughnessLayerResistance * roughnessLayerHeight * rVariables.previousRoughnessTemperature + timeStepSize *
            initialSoilTemperature + timeStepSize * currentWindSpeed * roughnessLayerResistance * surfaceRoughnessFactor *
            frictionDragCoefficient * frictionDragCoefficient * currentAirTemperature) / c;
        rVariables.roughnessTemperature += currentRoughnessTemperature / TNumNodes;
    }
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateNodalFluxes(
    const ProcessInfo& CurrentProcessInfo,
    ElementVariables& rVariables)
{
    const double albedoCoefficient = rVariables.albedoCoefficient;
    const double firstCoverStorageCoefficient = rVariables.firstCoverStorageCoefficient;
    const double secondCoverStorageCoefficient = rVariables.secondCoverStorageCoefficient;
    const double thirdCoverStorageCoefficient = rVariables.thirdCoverStorageCoefficient;
    const double buildEnvironmentRadiation = rVariables.buildEnvironmentRadiation;
    const double minimalStorage = rVariables.minimalStorage;
    const double maximalStorage = rVariables.maximalStorage;

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

    rVariables.previousStorage = rVariables.waterStorage;
    rVariables.previousRadiation = rVariables.netRadiation;
    rVariables.waterStorage = 0.0;
    rVariables.netRadiation = 0.0;

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
        const double sensibleHeatFluxRight = -airHeatCapacity * airDensity * rVariables.roughnessTemperature / roughnessLayerResistance;

        // Eq 5.35
        const double atmosphericResistance = 1.0 / (0.007 + 0.0056 * windSpeed);

        // Eq. 5.12
        const double saturatedVaporPressure = 6.11 * std::exp(17.27 * atmosphericTemperature / (atmosphericTemperature + 237.3));

        // Eq 5.13
        const double vaporPressureIncrement = 4098.0 * saturatedVaporPressure / (std::pow((atmosphericTemperature + 237.3), 2.0));

        // Eq 5.14
        const double actualVaporPressure = humidity / 100.0 * saturatedVaporPressure;

        // Eq 5.16
        const double shortWaveRadiation = (1.0 - albedoCoefficient) * incomingRadiation;

        // Eq 5.18
        const double emittedLongWaveRadiation = boltzmannCoefficient * std::pow(initialSoilTemperature + 273.15, 4.0);

        // Eq 5.17
        const double absorbedLongWaveRadiation = effectiveEmissivity * boltzmannCoefficient * std::pow(atmosphericTemperature + 273.15, 4.0);

        // Eq 5.15
        const double netRadiation = shortWaveRadiation + absorbedLongWaveRadiation - emittedLongWaveRadiation;

        // Eq 5.20
        const double surfaceHeatStorage = firstCoverStorageCoefficient * netRadiation + secondCoverStorageCoefficient * (netRadiation - rVariables.previousRadiation) /
            timeStepSize + thirdCoverStorageCoefficient;

        // Eq 5.34
        double latentHeatFlux = (vaporPressureIncrement * (netRadiation + buildEnvironmentRadiation - surfaceHeatStorage) + airHeatCapacity * airDensity *
            (saturatedVaporPressure - actualVaporPressure) / atmosphericResistance) / (vaporPressureIncrement + psychometricConstant *  // division is replaced to * based on (3.34)
                (1.0 + surfaceResistance / atmosphericResistance));   //Where is G (5.34)?
        latentHeatFlux = std::max(latentHeatFlux, 0.0);

        double actualEvaporation = 0.0;
        double actualPrecipitation = 0.0;

        // Eq 5.36
        const double potentialEvaporation = latentHeatFlux / (waterDensity * latentEvaporationHeat);
        const double potentialStorage = rVariables.previousStorage + timeStepSize * (precipitation - potentialEvaporation);
        if (potentialStorage > maximalStorage)
        {
            actualEvaporation = potentialEvaporation;
            actualPrecipitation = (maximalStorage - rVariables.previousStorage) / timeStepSize + actualEvaporation;
        }
        else if (potentialStorage < minimalStorage)
        {
            actualPrecipitation = precipitation;
            actualEvaporation = (rVariables.previousStorage - minimalStorage) / timeStepSize + actualPrecipitation;
        }
        else
        {
            actualEvaporation = potentialEvaporation;
            actualPrecipitation = precipitation;
        }
        const double actualStorage = rVariables.previousStorage + timeStepSize * (actualPrecipitation - actualEvaporation);
        latentHeatFlux = actualEvaporation * waterDensity * latentEvaporationHeat;

        // Eq 5.31
        double subsurfaceHeatFlux = netRadiation - sensibleHeatFluxRight - latentHeatFlux + buildEnvironmentRadiation - surfaceHeatStorage;

        rVariables.netRadiation += netRadiation / TNumNodes;
        rVariables.waterStorage += actualStorage / TNumNodes;
        rVariables.leftHandSideFlux[i] = sensibleHeatFluxLeft;
        rVariables.rightHandSideFlux[i] = subsurfaceHeatFlux;
    }
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;

    //Resetting the LHS
    if (rLeftHandSideMatrix.size1() != conditionSize)
        rLeftHandSideMatrix.resize(conditionSize, conditionSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(conditionSize, conditionSize);

    //Resetting the RHS
    if (rRightHandSideVector.size() != conditionSize)
        rRightHandSideVector.resize(conditionSize, false);
    noalias(rRightHandSideVector) = ZeroVector(conditionSize);

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::InitializeProperties()
{
    KRATOS_TRY

    const PropertiesType& rProp = this->GetProperties();

    rVariables.albedoCoefficient = rProp[ALPHA_COEFFICIENT];
    rVariables.firstCoverStorageCoefficient = rProp[A1_COEFFICIENT];
    rVariables.secondCoverStorageCoefficient = rProp[A2_COEFFICIENT];
    rVariables.thirdCoverStorageCoefficient = rProp[A3_COEFFICIENT];
    rVariables.buildEnvironmentRadiation = rProp[QF_COEFFICIENT];
    rVariables.minimalStorage = rProp[SMIN_COEFFICIENT];
    rVariables.maximalStorage = rProp[SMAX_COEFFICIENT];

    const GeometryType& Geom = this->GetGeometry();
    rVariables.roughnessTemperature = Geom[0].FastGetSolutionStepValue(AIR_TEMPERATURE, 1);   // This value is not read correctly
    rVariables.waterStorage = 0.0;                                                            // it is related to the initial value of the table
    rVariables.netRadiation = Geom[0].FastGetSolutionStepValue(SOLAR_RADIATION, 1);  // This value is not read correctly, initial value of the table

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TMicroClimateFluxCondition<TDim, TNumNodes>::CalculateKinematics(
    ElementVariables& rVariables,
    unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    rVariables.Np = row(rVariables.NContainer, PointNumber);
    rVariables.detJ = rVariables.detJContainer[PointNumber];

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
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
