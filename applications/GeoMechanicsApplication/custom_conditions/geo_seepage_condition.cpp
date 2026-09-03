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

#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/dof_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"
#include "includes/variables.h"

namespace Kratos
{

GeoSeepageCondition::GeoSeepageCondition() : GeoSeepageCondition(0, nullptr, nullptr) {}

GeoSeepageCondition::GeoSeepageCondition(IndexType ConditionId, GeometryType::Pointer pGeometry)
    : GeoSeepageCondition(ConditionId, std::move(pGeometry), nullptr)
{
}

GeoSeepageCondition::GeoSeepageCondition(IndexType               ConditionId,
                                         GeometryType::Pointer   pGeometry,
                                         PropertiesType::Pointer pProperties)
    : Condition(ConditionId, std::move(pGeometry), std::move(pProperties))
{
}

GeoSeepageCondition::~GeoSeepageCondition() = default;

Condition::Pointer GeoSeepageCondition::Create(IndexType               ConditionId,
                                               const NodesArrayType&   rNodes,
                                               PropertiesType::Pointer pProperties) const
{
    return Create(ConditionId, GetGeometry().Create(rNodes), pProperties);
}

Condition::Pointer GeoSeepageCondition::Create(IndexType               ConditionId,
                                               GeometryType::Pointer   pGeometry,
                                               PropertiesType::Pointer pProperties) const
{
    return make_intrusive<GeoSeepageCondition>(ConditionId, pGeometry, pProperties);
}

void GeoSeepageCondition::GetDofList(DofsVectorType& rResult, const ProcessInfo&) const
{
    rResult = GetDofs();
}

void GeoSeepageCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void GeoSeepageCondition::CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                                               Vector&            rRightHandSideVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void GeoSeepageCondition::CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo&)
{
    // A seepage condition never contributes to the system matrix: in Dirichlet mode the fixed
    // degrees of freedom are handled by the builder and solver, and in Neumann mode the flux is zero.
    const auto number_of_nodes = GetGeometry().PointsNumber();
    rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);
}

void GeoSeepageCondition::CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo&)
{
    const auto number_of_nodes = GetGeometry().PointsNumber();
    rRightHandSideVector.resize(number_of_nodes, false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);
}

Condition::DofsVectorType GeoSeepageCondition::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
}

int GeoSeepageCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto base_check_result = Condition::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF(GetGeometry().PointsNumber() < 2)
        << "GeoSeepageCondition " << Id() << " needs at least two nodes, but has "
        << GetGeometry().PointsNumber() << std::endl;

    for (const auto& r_node : GetGeometry()) {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(WATER_PRESSURE))
            << "Missing variable WATER_PRESSURE on node " << r_node.Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_node.HasDofFor(WATER_PRESSURE))
            << "Missing degree of freedom for WATER_PRESSURE on node " << r_node.Id() << std::endl;
    }

    return base_check_result;

    KRATOS_CATCH("")
}

std::string GeoSeepageCondition::Info() const { return "GeoSeepageCondition"; }

void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

} // namespace Kratos
