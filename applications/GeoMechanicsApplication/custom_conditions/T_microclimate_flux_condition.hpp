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
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

struct WaterFluxes
{
    double precipitation = 0.0;
    double evaporation = 0.0;
};

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TMicroClimateFluxCondition : public GeoTCondition<TDim,TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TMicroClimateFluxCondition);
    
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    TMicroClimateFluxCondition() : GeoTCondition<TDim, TNumNodes>() {}

    TMicroClimateFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry)
    {
    }

    TMicroClimateFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
        : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    Condition::Pointer Create(IndexType NewId,
                              const NodesArrayType& rNodes,
                              Properties::Pointer pProperties ) const override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

private:
    struct ElementVariables
    {
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        array_1d<double, TNumNodes> rightHandSideFlux;
    };

    void CalculateAndAddRHS(Vector& rRightHandSideVector,
                            const Vector& rNodalTemperatures, const array_1d<double, TNumNodes>& rLeftHandSideFluxes);
    void CalculateAndAddLHS(Matrix& rLeftHandSideMatrix, const array_1d<double, TNumNodes>& rLeftHandSideFluxes);

    double CalculateIntegrationCoefficient(const Matrix& rJacobian,
                                           double Weight);

    void CalculateRoughness(const ProcessInfo& rCurrentProcessInfo);

    void InitializeProperties();

    double CalculateNetRadiation(unsigned int index);

    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }

    bool mIsInitialised = false;
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
    ElementVariables mVariables;

    array_1d<double, TNumNodes> CalculateLeftHandSideFluxes();
    double CalculatePotentialEvaporation(unsigned int i,
                                   const double net_radiation,
                                   const double surface_heat_storage);
    void SetWaterStorage(double time_step_size,
                          double previous_storage,
                          double actual_precipitation,
                          double actual_evaporation);
    void SetRightHandSideFlux(unsigned int i,
                              const double net_radiation,
                              const double surface_heat_storage,
                              double actual_evaporation);
    void SetNetRadiation();
    WaterFluxes CalculateWaterFluxes(unsigned int i,
                                     const double time_step_size,
                                     const double previous_storage,
                                     const double net_radiation,
                                     const double surface_heat_storage);
    double CalculateSurfaceHeatStorage(const double time_step_size,
                                       const double previous_radiation,
                                       const double net_radiation) const;
    void SetWaterStorage(const double time_step_size,
                         const double previous_storage,
                         const double previous_radiation);
    void SetRightHandSideFluxes(const double time_step_size,
                                const double previous_storage,
                                const double previous_radiation);
}; // class TMicroClimateFluxCondition.

} // namespace Kratos.