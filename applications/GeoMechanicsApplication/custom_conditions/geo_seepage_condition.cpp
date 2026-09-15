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
#include "includes/serializer.h"
#include "includes/variables.h"

using namespace std::string_literals;

namespace Kratos
{

GeoSeepageCondition::GeoSeepageCondition(IndexType               ConditionId,
                                         GeometryType::Pointer   pGeometry,
                                         PropertiesType::Pointer pProperties)
    : Condition(ConditionId, std::move(pGeometry), std::move(pProperties))
{
}

Condition::Pointer GeoSeepageCondition::Create(IndexType               ConditionId,
                                               const NodesArrayType&   rNodes,
                                               PropertiesType::Pointer pProperties) const
{
    return Create(ConditionId, GetGeometry().Create(rNodes), std::move(pProperties));
}

Condition::Pointer GeoSeepageCondition::Create(IndexType               ConditionId,
                                               GeometryType::Pointer   pGeometry,
                                               PropertiesType::Pointer pProperties) const
{
    return make_intrusive<GeoSeepageCondition>(ConditionId, pGeometry, pProperties);
}

int GeoSeepageCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
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
}

std::string GeoSeepageCondition::Info() const { return "GeoSeepageCondition"s; }

void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

} // namespace Kratos
