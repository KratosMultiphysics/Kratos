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
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    // Default constructor
    TMicroClimateFluxCondition() : GeoTCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    TMicroClimateFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : GeoTCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    TMicroClimateFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : GeoTCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

private:
    struct ElementVariables
    {
        double albedoCoefficient;
        double firstCoverStorageCoefficient;
        double secondCoverStorageCoefficient;
        double thirdCoverStorageCoefficient;
        double buildEnvironmentRadiation;
        double minimalStorage;
        double maximalStorage;

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
    ElementVariables mVariables;
}; // class TMicroClimateFluxCondition.

} // namespace Kratos.