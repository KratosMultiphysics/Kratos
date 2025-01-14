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

// Project includes
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include <custom_utilities/variables_utilities.hpp>

#include "includes/variables.h"

namespace Kratos
{

// Default Constructor
LineNormalLoad2DDiffOrderCondition::LineNormalLoad2DDiffOrderCondition()
    : LineLoad2DDiffOrderCondition()
{
}

// Constructor 1
LineNormalLoad2DDiffOrderCondition::LineNormalLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : LineLoad2DDiffOrderCondition(NewId, pGeometry)
{
}

// Constructor 2
LineNormalLoad2DDiffOrderCondition::LineNormalLoad2DDiffOrderCondition(IndexType NewId,
                                                                       GeometryType::Pointer pGeometry,
                                                                       PropertiesType::Pointer pProperties)
    : LineLoad2DDiffOrderCondition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer LineNormalLoad2DDiffOrderCondition::Create(IndexType             NewId,
                                                              NodesArrayType const& ThisNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LineNormalLoad2DDiffOrderCondition::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return make_intrusive<LineNormalLoad2DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void LineNormalLoad2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                  unsigned int        PointNumber)
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();

    Vector tangential_vector = ZeroVector(3);
    tangential_vector[0]     = column(rVariables.JContainer[PointNumber], 0)[0];
    tangential_vector[1]     = column(rVariables.JContainer[PointNumber], 0)[1];
    Vector tangential_stresses(r_geometry.PointsNumber());
    VariablesUtilities::GetNodalValues(r_geometry, TANGENTIAL_CONTACT_STRESS, tangential_stresses.begin());
    auto tangential_stress = MathUtils<>::Dot(rVariables.Nu, tangential_stresses);

    Vector out_of_plane_vector = ZeroVector(3);
    out_of_plane_vector[2]     = 1.0;
    Vector normal_vector       = ZeroVector(3);
    MathUtils<double>::CrossProduct(normal_vector, out_of_plane_vector, tangential_vector);
    Vector normal_stresses(r_geometry.PointsNumber());
    VariablesUtilities::GetNodalValues(r_geometry, NORMAL_CONTACT_STRESS, normal_stresses.begin());
    auto normal_stress = MathUtils<>::Dot(rVariables.Nu, normal_stresses);

    auto traction_vector = tangential_stress * tangential_vector + normal_stress * normal_vector;

    rVariables.ConditionVector.resize(2, false);
    rVariables.ConditionVector[0] = traction_vector[0];
    rVariables.ConditionVector[1] = traction_vector[1];

    KRATOS_CATCH("")
}

double LineNormalLoad2DDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    return IntegrationPoints[PointNumber].Weight();
}

void LineNormalLoad2DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                       ConditionVariables& rVariables)
{
    for (SizeType node = 0; node < this->GetGeometry().PointsNumber(); ++node) {
        rRightHandSideVector[2 * node] +=
            rVariables.Nu[node] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[2 * node + 1] +=
            rVariables.Nu[node] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }
}

std::string LineNormalLoad2DDiffOrderCondition::Info() const
{
    return "LineNormalLoad2DDiffOrderCondition";
}

} // Namespace Kratos.
