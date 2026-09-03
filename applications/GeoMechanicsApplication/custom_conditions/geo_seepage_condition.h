// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Wijtze Pieter Kikstra

#pragma once

#include <string>

#include "includes/condition.h"
#include "includes/define.h"

namespace Kratos
{

class Serializer;

// A seepage boundary condition on the WATER_PRESSURE degree of freedom.
//
// This condition holds no state of its own. A node's WATER_PRESSURE fixity is the single source of
// truth for its boundary type: fixed means a Dirichlet boundary at zero pressure, free means a
// zero-flux Neumann boundary. GeoSeepageNewtonRaphsonStrategy switches individual nodes while
// iterating; this condition only marks which element boundaries belong to the seepage face.
//
// The condition never contributes to the linear system.
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoSeepageCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoSeepageCondition);

    GeoSeepageCondition();
    GeoSeepageCondition(IndexType ConditionId, GeometryType::Pointer pGeometry);
    GeoSeepageCondition(IndexType ConditionId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    ~GeoSeepageCondition() override;

    Condition::Pointer Create(IndexType               ConditionId,
                              const NodesArrayType&   rNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType               ConditionId,
                              GeometryType::Pointer   pGeometry,
                              PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rResult, const ProcessInfo&) const override;
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo&) override;
    void CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo&) override;

    [[nodiscard]] int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    friend Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
