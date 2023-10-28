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
//                   John van Esch
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/T_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TNormalFluxCondition : public TCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TNormalFluxCondition);
    
    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    
    // ============================================================================================
    // ============================================================================================

    // Default constructor
    TNormalFluxCondition();
    
    // Constructor 1
    TNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    
    // Constructor 2
    TNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    ~TNormalFluxCondition() override;

    // ============================================================================================
    // ============================================================================================

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
    // ============================================================================================
    // ============================================================================================

protected:

    struct NormalFluxVariables
    {
        double normalFlux;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> N;
        array_1d<double,TNumNodes> fluxVector;
    };
    
    // ============================================================================================
    // ============================================================================================
                                    
    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables);

    virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient,
        const Matrix& Jacobian,
        const double& Weight);
    
    // ============================================================================================
    // ============================================================================================

private:

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
    
}; // class TNormalFluxCondition.

} // namespace Kratos.