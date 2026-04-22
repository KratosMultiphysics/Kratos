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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwNormalFaceLoadCondition : public UPwCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwNormalFaceLoadCondition);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    UPwNormalFaceLoadCondition();

    UPwNormalFaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    UPwNormalFaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    std::string Info() const override;

protected:
    struct NormalFaceLoadVariables {
        array_1d<double, TNumNodes> NormalStressVector;
        array_1d<double, TNumNodes> TangentialStressVector;
    };

    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

    void InitializeConditionVariables(NormalFaceLoadVariables& rVariables, const GeometryType& Geom);

    void CalculateTractionVector(array_1d<double, TDim>&        rTractionVector,
                                 const Matrix&                  Jacobian,
                                 const Matrix&                  NContainer,
                                 const NormalFaceLoadVariables& Variables,
                                 const unsigned int&            GPoint);

    virtual double CalculateIntegrationCoefficient(const IndexType PointNumber,
                                                   const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class UPwNormalFaceLoadCondition.

} // namespace Kratos.