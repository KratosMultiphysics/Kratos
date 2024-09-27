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
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

// Default Constructor
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition()
    : SurfaceLoad3DDiffOrderCondition()
{
}

//----------------------------------------------------------------------------------------
// Constructor 1
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId,
                                                                             GeometryType::Pointer pGeometry)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
// Constructor 2
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId,
                                                                             GeometryType::Pointer pGeometry,
                                                                             PropertiesType::Pointer pProperties)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Condition::Pointer SurfaceNormalLoad3DDiffOrderCondition::Create(IndexType             NewId,
                                                                 NodesArrayType const& ThisNodes,
                                                                 PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalLoad3DDiffOrderCondition(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void SurfaceNormalLoad3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                     unsigned int PointNumber)
{
    KRATOS_TRY

    Vector normal_vector(3);
    MathUtils<double>::CrossProduct(normal_vector, column(rVariables.JContainer[PointNumber], 0),
                                    column(rVariables.JContainer[PointNumber], 1));

    const auto& r_geometry = GetGeometry();

    Vector normal_stresses(r_geometry.PointsNumber());
    std::transform(r_geometry.begin(), r_geometry.end(), normal_stresses.begin(), [](const auto& r_node) {
        return r_node.FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    });

    // Since the normal vector is pointing outwards for the 3D conditions, the normal stress
    // should switch sign, such that positive normal contact stress is defined inwards.
    const double normal_stress = -1 * MathUtils<>::Dot(rVariables.Nu, normal_stresses);
    rVariables.ConditionVector = normal_stress * normal_vector;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
double SurfaceNormalLoad3DDiffOrderCondition::CalculateIntegrationCoefficient(
    const IndexType                                 PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    return IntegrationPoints[PointNumber].Weight();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void SurfaceNormalLoad3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                                                          ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    SizeType       Index;

    for (SizeType i = 0; i < NumUNodes; ++i) {
        Index = i * 3;

        rRightHandSideVector[Index] +=
            rVariables.Nu[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index + 1] +=
            rVariables.Nu[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index + 2] +=
            rVariables.Nu[i] * rVariables.ConditionVector[2] * rVariables.IntegrationCoefficient;
    }
}

std::string SurfaceNormalLoad3DDiffOrderCondition::Info() const
{
    return "SurfaceNormalLoad3DDiffOrderCondition";
}

} // Namespace Kratos.
