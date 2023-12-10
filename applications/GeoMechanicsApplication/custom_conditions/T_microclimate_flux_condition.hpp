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
        double roughnessTemperature = 0.0;
        double netRadiation = 0.0;
        double waterStorage = 0.0;

        array_1d<double, TNumNodes> leftHandSideFlux;
        array_1d<double, TNumNodes> rightHandSideFlux;
    };

    void CalculateAll(Matrix&            rLeftHandSideMatrix,
                      Vector&            rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateAndAddRHS(Vector& rRightHandSideVector,
                            const Vector& rNodalTemperatures);
    void CalculateAndAddLHS(Matrix& rLeftHandSideMatrix);

    double CalculateIntegrationCoefficient(const Matrix& Jacobian,
                                           double Weight);

    void CalculateRoughness(const ProcessInfo& CurrentProcessInfo);

    void CalculateNodalFluxes(const ProcessInfo& CurrentProcessInfo);

    void InitializeProperties();

    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }

    bool mIsInitialised = false;
    double mAlbedoCoefficient = 0.0;
    double mFirstCoverStorageCoefficient = 0.0;
    double mSecondCoverStorageCoefficient = 0.0;
    double mThirdCoverStorageCoefficient = 0.0;
    double mBuildEnvironmentRadiation = 0.0;
    double mMinimalStorage = 0.0;
    double mMaximalStorage = 0.0;
    ElementVariables mVariables;
}; // class TMicroClimateFluxCondition.

} // namespace Kratos.