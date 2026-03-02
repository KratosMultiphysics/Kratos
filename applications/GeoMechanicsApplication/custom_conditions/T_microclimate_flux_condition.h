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
//  Main authors:    Mohamed Nabi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/T_condition.h"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

struct WaterFluxes {
    double precipitation = 0.0;
    double evaporation   = 0.0;
};

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTMicroClimateFluxCondition
    : public GeoTCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTMicroClimateFluxCondition);

    using IndexType      = std::size_t;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    GeoTMicroClimateFluxCondition();

    GeoTMicroClimateFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    GeoTMicroClimateFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    Condition::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, Properties::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override;

private:
    void InitializeProperties();

    void CalculateAndAddLHS(Matrix&                            rLeftHandSideMatrix,
                            const array_1d<double, TNumNodes>& rN,
                            double                             IntegrationCoefficient,
                            const array_1d<double, TNumNodes>& rLeftHandSideFluxes);

    void CalculateAndAddRHS(Vector&                            rRightHandSideVector,
                            const array_1d<double, TNumNodes>& rN,
                            double                             IntegrationCoefficient,
                            const Vector&                      rNodalTemperatures,
                            const array_1d<double, TNumNodes>& rLeftHandSideFluxes,
                            const array_1d<double, TNumNodes>& rRightHandSideFluxes);

    array_1d<double, TNumNodes> CalculateLeftHandSideFluxes() const;
    array_1d<double, TNumNodes> CalculateRightHandSideFluxes(double TimeStepSize,
                                                             double PreviousStorage,
                                                             double PreviousRadiation) const;
    double CalculateRightHandSideFlux(double NetRadiation, double SurfaceHeatStorage, double ActualEvaporation) const;

    double CalculateCurrentWaterStorage(double TimeStepSize, double PreviousStorage, double PreviousRadiation) const;

    double CalculateCurrentNetRadiation() const;
    double CalculateNetRadiation(unsigned int NodeIndex) const;

    double CalculateSurfaceHeatStorage(double TimeStepSize, double PreviousRadiation, double NetRadiation) const;

    WaterFluxes CalculateWaterFluxes(unsigned int NodeIndex,
                                     double       TimeStepSize,
                                     double       PreviousStorage,
                                     double       NetRadiation,
                                     double       SurfaceHeatStorage) const;

    double CalculatePotentialEvaporation(unsigned int NodeIndex, double NetRadiation, double SurfaceHeatStorage) const;

    void CalculateRoughness(const ProcessInfo& rCurrentProcessInfo);

    double CalculateSurfaceRoughnessFactor(double CurrentAirTemperature,
                                           double PreviousRoughnessTemperature,
                                           double RichardsonBulkModulus,
                                           double FrictionDragCoefficient) const;

    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    bool   mIsInitialized                 = false;
    double mAlbedoCoefficient             = 0.0;
    double mFirstCoverStorageCoefficient  = 0.0;
    double mSecondCoverStorageCoefficient = 0.0;
    double mThirdCoverStorageCoefficient  = 0.0;
    double mBuildEnvironmentRadiation     = 0.0;
    double mMinimalStorage                = 0.0;
    double mMaximalStorage                = 0.0;
    double mRoughnessTemperature          = 0.0;
    double mNetRadiation                  = 0.0;
    double mWaterStorage                  = 0.0;
    double mWaterDensity                  = 0.0;
};

// class TMicroClimateFluxCondition.

} // namespace Kratos.