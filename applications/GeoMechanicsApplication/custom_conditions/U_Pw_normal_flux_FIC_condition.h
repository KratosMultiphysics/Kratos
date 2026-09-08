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
#include "custom_conditions/U_Pw_normal_flux_condition.h"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwNormalFluxFICCondition
    : public UPwNormalFluxCondition<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwNormalFluxFICCondition);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using GeometryType   = Geometry<Node>;
    using NodesArrayType = GeometryType::PointsArrayType;

    UPwNormalFluxFICCondition();

    UPwNormalFluxFICCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    UPwNormalFluxFICCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    std::string Info() const override;

protected:
    struct NormalFluxVariables {
        double                      IntegrationCoefficient;
        array_1d<double, TNumNodes> Np;
    };

    struct NormalFluxFICVariables {
        double DtPressureCoefficient;
        double ElementLength;
        double BiotModulusInverse;

        array_1d<double, TNumNodes> DtPressureVector;
    };

    void CalculateAll(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);

    void CalculateAndAddLHSStabilization(Matrix&                 rLeftHandSideMatrix,
                                         NormalFluxVariables&    rVariables,
                                         NormalFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassMatrix(Matrix&                       rLeftHandSideMatrix,
                                           const NormalFluxVariables&    rVariables,
                                           const NormalFluxFICVariables& rFICVariables);

    void CalculateAndAddRHSStabilization(Vector&                 rRightHandSideVector,
                                         NormalFluxVariables&    rVariables,
                                         NormalFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassFlow(Vector&                       rRightHandSideVector,
                                         NormalFluxVariables&          rVariables,
                                         const NormalFluxFICVariables& rFICVariables);

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

}; // class UPwNormalFluxFICCondition.

} // namespace Kratos.
