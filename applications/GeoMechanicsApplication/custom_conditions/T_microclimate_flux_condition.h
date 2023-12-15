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
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

struct WaterFluxes
{
    double precipitation = 0.0;
    double evaporation = 0.0;
};

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTMicroClimateFluxCondition
    : public GeoTCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTMicroClimateFluxCondition);

    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    GeoTMicroClimateFluxCondition() : GeoTCondition<TDim, TNumNodes>() {}

    GeoTMicroClimateFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry)
    {
    }

    GeoTMicroClimateFluxCondition(IndexType NewId,
                                  GeometryType::Pointer pGeometry,
                                  Properties::Pointer pProperties)
        : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    Condition::Pointer Create(IndexType NewId,
                              const NodesArrayType& rNodes,
                              Properties::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

private:
    void InitializeProperties();

    void CalculateAndAddLHS(Matrix& rLeftHandSideMatrix,
                            const array_1d<double, TNumNodes>& rN,
                            double IntegrationCoefficient,
                            const array_1d<double, TNumNodes>& rLeftHandSideFluxes);

    void CalculateAndAddRHS(Vector& rRightHandSideVector,
                            const array_1d<double, TNumNodes>& rN,
                            double IntegrationCoefficient,
                            const Vector& rNodalTemperatures,
                            const array_1d<double, TNumNodes>& rLeftHandSideFluxes,
                            const array_1d<double, TNumNodes>& rRightHandSideFluxes);

    array_1d<double, TNumNodes> CalculateLeftHandSideFluxes() const;
    array_1d<double, TNumNodes> CalculateRightHandSideFluxes(double TimeStepSize,
                                                             double PreviousStorage,
                                                             double PreviousRadiation) const;
    double CalculateRightHandSideFlux(double NetRadiation,
                                      double SurfaceHeatStorage,
                                      double ActualEvaporation) const;

    double CalculateCurrentWaterStorage(double TimeStepSize,
                                        double PreviousStorage,
                                        double PreviousRadiation) const;

    double CalculateCurrentNetRadiation() const;
    double CalculateNetRadiation(unsigned int NodeIndex) const;

    double CalculateSurfaceHeatStorage(double TimeStepSize,
                                       double PreviousRadiation,
                                       double NetRadiation) const;

    WaterFluxes CalculateWaterFluxes(unsigned int NodeIndex,
                                     double TimeStepSize,
                                     double PreviousStorage,
                                     double NetRadiation,
                                     double SurfaceHeatStorage) const;

    double CalculatePotentialEvaporation(unsigned int NodeIndex,
                                         double NetRadiation,
                                         double SurfaceHeatStorage) const;

    void CalculateRoughness(const ProcessInfo& rCurrentProcessInfo);

    double CalculateSurfaceRoughnessFactor(double CurrentAirTemperature,
                                           double PreviousRoughnessTemperature,
                                           double RichardsonBulkModulus,
                                           double FrictionDragCoefficient) const;

    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
        rSerializer.save("mIsInitialized", mIsInitialized);
        rSerializer.save("mAlbedoCoefficient", mAlbedoCoefficient);
        rSerializer.save("mFirstCoverStorageCoefficient", mFirstCoverStorageCoefficient);
        rSerializer.save("mSecondCoverStorageCoefficient", mSecondCoverStorageCoefficient);
        rSerializer.save("mThirdCoverStorageCoefficient", mThirdCoverStorageCoefficient);
        rSerializer.save("mBuildEnvironmentRadiation", mBuildEnvironmentRadiation);
        rSerializer.save("mMinimalStorage", mMinimalStorage);
        rSerializer.save("mMaximalStorage", mMaximalStorage);
        rSerializer.save("mRoughnessTemperature", mRoughnessTemperature);
        rSerializer.save("mNetRadiation", mNetRadiation);
        rSerializer.save("mWaterStorage", mWaterStorage);
        rSerializer.save("mWaterDensity", mWaterDensity);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
        rSerializer.load("mIsInitialized", mIsInitialized);
        rSerializer.load("mAlbedoCoefficient", mAlbedoCoefficient);
        rSerializer.load("mFirstCoverStorageCoefficient", mFirstCoverStorageCoefficient);
        rSerializer.load("mSecondCoverStorageCoefficient", mSecondCoverStorageCoefficient);
        rSerializer.load("mThirdCoverStorageCoefficient", mThirdCoverStorageCoefficient);
        rSerializer.load("mBuildEnvironmentRadiation", mBuildEnvironmentRadiation);
        rSerializer.load("mMinimalStorage", mMinimalStorage);
        rSerializer.load("mMaximalStorage", mMaximalStorage);
        rSerializer.load("mRoughnessTemperature", mRoughnessTemperature);
        rSerializer.load("mNetRadiation", mNetRadiation);
        rSerializer.load("mWaterStorage", mWaterStorage);
        rSerializer.load("mWaterDensity", mWaterDensity);
    }

    bool mIsInitialized = false;
    double mAlbedoCoefficient = 0.0;
    double mFirstCoverStorageCoefficient = 0.0;
    double mSecondCoverStorageCoefficient = 0.0;
    double mThirdCoverStorageCoefficient = 0.0;
    double mBuildEnvironmentRadiation = 0.0;
    double mMinimalStorage = 0.0;
    double mMaximalStorage = 0.0;
    double mRoughnessTemperature = 0.0;
    double mNetRadiation = 0.0;
    double mWaterStorage = 0.0;
    double mWaterDensity = 0.0;
};

// class TMicroClimateFluxCondition.

} // namespace Kratos.