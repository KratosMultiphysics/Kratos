// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/dof_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/variables.h"

namespace Kratos
{

GeoSeepageCondition::GeoSeepageCondition() : GeoSeepageCondition(0, nullptr, nullptr) {}

GeoSeepageCondition::GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoSeepageCondition(NewId, pGeometry, nullptr)
{
}

GeoSeepageCondition::GeoSeepageCondition(IndexType               NewId,
                                         GeometryType::Pointer   pGeometry,
                                         PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

GeoSeepageCondition::~GeoSeepageCondition() = default;

Condition::Pointer GeoSeepageCondition::Create(IndexType               NewId,
                                               NodesArrayType const&   rThisNodes,
                                               PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

Condition::Pointer GeoSeepageCondition::Create(IndexType               NewId,
                                               GeometryType::Pointer   pGeom,
                                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoSeepageCondition>(NewId, pGeom, pProperties);
}

void GeoSeepageCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

void GeoSeepageCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void GeoSeepageCondition::CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                                               Vector&            rRightHandSideVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void GeoSeepageCondition::CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo&)
{
    KRATOS_TRY

    // A seepage condition never contributes to the system matrix: in Dirichlet mode the fixed
    // degrees of freedom are handled by the builder and solver, and in Neumann mode the flux is zero.
    const auto number_of_nodes = GetGeometry().PointsNumber();
    if (rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

    KRATOS_CATCH("")
}

void GeoSeepageCondition::CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo&)
{
    KRATOS_TRY

    const auto number_of_nodes = GetGeometry().PointsNumber();
    if (rRightHandSideVector.size() != number_of_nodes) {
        rRightHandSideVector.resize(number_of_nodes, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);

    KRATOS_CATCH("")
}

Condition::DofsVectorType GeoSeepageCondition::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
}

void GeoSeepageCondition::Initialize(const ProcessInfo&)
{
    KRATOS_TRY

    // Properties are normally shared by many conditions. Since the boundary type is switched per
    // condition, each condition needs its own copy. The guard keeps this idempotent, so that a
    // second Initialize cannot undo a switch that has already been made.
    if (mHasOwnProperties) return;

    SetProperties(Kratos::make_shared<PropertiesType>(GetProperties()));
    mHasOwnProperties = true;

    KRATOS_CATCH("")
}

void GeoSeepageCondition::InitializeNonLinearIteration(const ProcessInfo&)
{
    KRATOS_TRY

    const auto& r_boundary_type = GetBoundaryType();
    KRATOS_ERROR_IF(r_boundary_type != Geo::SeepageDirichletType && r_boundary_type != Geo::SeepageNeumannType)
        << "Unknown seepage boundary type '" << r_boundary_type << "' for GeoSeepageCondition "
        << Id() << ", expected '" << Geo::SeepageDirichletType << "' or '"
        << Geo::SeepageNeumannType << "'" << std::endl;

    const auto is_dirichlet = r_boundary_type == Geo::SeepageDirichletType;
    for (auto& r_node : GetGeometry()) {
        if (is_dirichlet) {
            // A seepage face lets water out at atmospheric pressure.
            r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
            r_node.Fix(WATER_PRESSURE);
        } else {
            // Zero flux: freeing the degree of freedom is all that is needed, since a zero flux
            // makes no contribution to the system.
            r_node.Free(WATER_PRESSURE);
        }
    }

    KRATOS_CATCH("")
}

const std::string& GeoSeepageCondition::GetBoundaryType() const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(GEO_SEEPAGE_BOUNDARY_TYPE))
        << "GEO_SEEPAGE_BOUNDARY_TYPE is not defined for GeoSeepageCondition " << Id() << std::endl;

    return GetProperties()[GEO_SEEPAGE_BOUNDARY_TYPE];
}

void GeoSeepageCondition::SetBoundaryType(const std::string& rBoundaryType)
{
    KRATOS_ERROR_IF(rBoundaryType != Geo::SeepageDirichletType && rBoundaryType != Geo::SeepageNeumannType)
        << "Unknown seepage boundary type '" << rBoundaryType << "' for GeoSeepageCondition " << Id()
        << ", expected '" << Geo::SeepageDirichletType << "' or '" << Geo::SeepageNeumannType << "'"
        << std::endl;

    GetProperties().SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, rBoundaryType);
}

std::string GeoSeepageCondition::Info() const { return "GeoSeepageCondition"; }

void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    rSerializer.save("HasOwnProperties", mHasOwnProperties);
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    rSerializer.load("HasOwnProperties", mHasOwnProperties);
}

} // namespace Kratos





