// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "includes/element.h"
#include "stress_state_policy.h"
#include "integration_scheme.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LineInterfaceElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LineInterfaceElement);

    using Element::GeometryType;
    using Element::PropertiesType;

    LineInterfaceElement();

    LineInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

private:
    Element::DofsVectorType GetDofs() const;

    std::unique_ptr<IntegrationScheme> mIntegrationScheme;
    std::unique_ptr<StressStatePolicy> mStressStatePolicy;
};

} // namespace Kratos
