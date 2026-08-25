// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <string>

#include "includes/condition.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

// A seepage boundary condition on the WATER_PRESSURE degree of freedom.
//
// This condition holds no mode of its own. A node's WATER_PRESSURE fixity is the single source of
// truth for its boundary type: fixed means a Dirichlet boundary at zero pressure, free means a
// zero-flux Neumann boundary. GeoSeepageNewtonRaphsonStrategy switches individual nodes while
// iterating; this condition only marks which nodes belong to the seepage face and gives them their
// initial state.
//
// The condition never contributes to the linear system.
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoSeepageCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoSeepageCondition);

    GeoSeepageCondition();
    GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    ~GeoSeepageCondition() override;

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   rThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType               NewId,
                              GeometryType::Pointer   pGeom,
                              PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    // Puts the seepage face into its initial, over-prescribed state: every node is a Dirichlet
    // boundary at zero pressure. The strategy releases nodes from there, one per iteration.
    void InitializeSolutionStep(const ProcessInfo&) override;

    [[nodiscard]] int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos








