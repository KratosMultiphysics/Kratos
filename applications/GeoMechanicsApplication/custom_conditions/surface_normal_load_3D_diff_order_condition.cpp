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
#include "custom_conditions/surface_normal_load_3D_diff_order_condition.hpp"
#include <custom_utilities/variables_utilities.hpp>

#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

// Default Constructor
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition()
    : SurfaceLoad3DDiffOrderCondition()
{
}

// Constructor 1
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId,
                                                                             GeometryType::Pointer pGeometry)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry)
{
}

// Constructor 2
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId,
                                                                             GeometryType::Pointer pGeometry,
                                                                             PropertiesType::Pointer pProperties)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer SurfaceNormalLoad3DDiffOrderCondition::Create(IndexType             NewId,
                                                                 NodesArrayType const& ThisNodes,
                                                                 PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer SurfaceNormalLoad3DDiffOrderCondition::Create(IndexType             NewId,
                                                                 GeometryType::Pointer pGeom,
                                                                 PropertiesType::Pointer pProperties) const
{
    return make_intrusive<SurfaceNormalLoad3DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void SurfaceNormalLoad3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                     unsigned int PointNumber)
{
    KRATOS_TRY

    Vector normal_vector(3);
    MathUtils<double>::CrossProduct(normal_vector, column(rVariables.JContainer[PointNumber], 0),
                                    column(rVariables.JContainer[PointNumber], 1));

    const auto& r_geometry = GetGeometry();
    Vector      normal_stresses(r_geometry.PointsNumber());
    VariablesUtilities::GetNodalValues(r_geometry, NORMAL_CONTACT_STRESS, normal_stresses.begin());

    // Since the normal vector is pointing outwards for the 3D conditions, the normal stress
    // should switch sign, such that positive normal contact stress is defined inwards.
    const auto normal_stress   = -1.0 * MathUtils<>::Dot(rVariables.Nu, normal_stresses);
    rVariables.ConditionVector = normal_stress * normal_vector;

    KRATOS_CATCH("")
}

double SurfaceNormalLoad3DDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              rJContainer,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints) const
{
    KRATOS_TRY

    return rIntegrationPoints[PointNumber].Weight();

    KRATOS_CATCH("")
}

void SurfaceNormalLoad3DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                          ConditionVariables& rVariables)
{
    for (SizeType node = 0; node < GetGeometry().PointsNumber(); ++node) {
        rRightHandSideVector[3 * node] +=
            rVariables.Nu[node] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[3 * node + 1] +=
            rVariables.Nu[node] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[3 * node + 2] +=
            rVariables.Nu[node] * rVariables.ConditionVector[2] * rVariables.IntegrationCoefficient;
    }
}

std::string SurfaceNormalLoad3DDiffOrderCondition::Info() const
{
    return "SurfaceNormalLoad3DDiffOrderCondition";
}

} // Namespace Kratos.
