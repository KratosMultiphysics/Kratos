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

namespace Kratos::Geo
{

// The only legal values of GEO_SEEPAGE_BOUNDARY_TYPE.
inline const std::string SeepageDirichletType = "Dirichlet";
inline const std::string SeepageNeumannType   = "Neumann";

} // namespace Kratos::Geo

namespace Kratos
{

// A mixed seepage boundary condition on the WATER_PRESSURE degree of freedom.
//
// The condition acts either as a Dirichlet boundary (WATER_PRESSURE fixed at zero) or as a
// zero-flux Neumann boundary (WATER_PRESSURE free). The mode is stored in the condition's
// Properties under GEO_SEEPAGE_BOUNDARY_TYPE and is re-read on every non-linear iteration, so
// it can be switched from outside the condition while solving.
//
// The condition never contributes to the linear system; it only changes the fixity of its nodes.
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

    // Gives this condition its own copy of its Properties, so that its boundary type can be
    // switched independently of the other conditions that originally shared those Properties.
    void Initialize(const ProcessInfo&) override;

    // Applies the currently configured boundary type to the nodes of this condition. Called once
    // per non-linear iteration, so that a boundary type switched from outside takes effect here.
    void InitializeNonLinearIteration(const ProcessInfo&) override;

    // Returns the currently configured boundary type. Errors if GEO_SEEPAGE_BOUNDARY_TYPE is absent.
    [[nodiscard]] const std::string& GetBoundaryType() const;

    // Sets the boundary type on this condition's Properties. Errors on any value other than
    // Geo::SeepageDirichletType or Geo::SeepageNeumannType.
    void SetBoundaryType(const std::string& rBoundaryType);

    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    bool mHasOwnProperties = false;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos





