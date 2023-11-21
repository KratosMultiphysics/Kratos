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

#include "custom_conditions/T_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/condition.h"
#include "includes/serializer.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TNormalFluxCondition
    : public TCondition<TDim, TNumNodes> {
public:
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using PropertiesType = Properties;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TNormalFluxCondition);

    TNormalFluxCondition();

    TNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    TNormalFluxCondition(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties);

    ~TNormalFluxCondition() override;

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& rThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new TNormalFluxCondition(
            NewId, this->GetGeometry().Create(rThisNodes), pProperties));
    }

protected:
    struct NormalFluxVariables {
        double NormalFlux;
        double IntegrationCoefficient;
        array_1d<double, TNumNodes> N;
        array_1d<double, TNumNodes> FluxVector;
    };

    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables);

    virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient,
                                                 const Matrix& Jacobian,
                                                 double Weight);

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    }
};

} // namespace Kratos
