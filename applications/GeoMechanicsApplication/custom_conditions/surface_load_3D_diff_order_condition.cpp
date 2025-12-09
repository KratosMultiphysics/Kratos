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
#include "custom_conditions/surface_load_3D_diff_order_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include <custom_utilities/variables_utilities.hpp>

namespace Kratos
{

SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition() : GeneralUPwDiffOrderCondition()
{
}

SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeneralUPwDiffOrderCondition(NewId, std::move(pGeometry))
{
}

SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition(IndexType             NewId,
                                                                 GeometryType::Pointer pGeometry,
                                                                 PropertiesType::Pointer pProperties)
    : GeneralUPwDiffOrderCondition(NewId, std::move(pGeometry), std::move(pProperties))
{
}

Condition::Pointer SurfaceLoad3DDiffOrderCondition::Create(IndexType             NewId,
                                                           NodesArrayType const& ThisNodes,
                                                           PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer SurfaceLoad3DDiffOrderCondition::Create(IndexType             NewId,
                                                           GeometryType::Pointer pGeom,
                                                           PropertiesType::Pointer pProperties) const
{
    return make_intrusive<SurfaceLoad3DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void SurfaceLoad3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const GeometryType& r_geom = GetGeometry();
    rVariables.ConditionVector.resize(3, false);
    noalias(rVariables.ConditionVector) = ZeroVector(3);

    for (SizeType node = 0; node < r_geom.PointsNumber(); ++node) {
        rVariables.ConditionVector += rVariables.Nu[node] * r_geom[node].FastGetSolutionStepValue(SURFACE_LOAD);
    }

    KRATOS_CATCH("")
}

double SurfaceLoad3DDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    KRATOS_TRY
    return ConditionUtilities::CalculateIntegrationCoefficient(
        JContainer[PointNumber], IntegrationPoints[PointNumber].Weight());
    KRATOS_CATCH("")
}

void SurfaceLoad3DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
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

void SurfaceLoad3DDiffOrderCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeneralUPwDiffOrderCondition)
}

void SurfaceLoad3DDiffOrderCondition::load(Serializer& rSerializer){KRATOS_SERIALIZE_LOAD_BASE_CLASS(
    rSerializer, GeneralUPwDiffOrderCondition)} std::string SurfaceLoad3DDiffOrderCondition::Info() const
{
    return "SurfaceLoad3DDiffOrderCondition";
}

} // Namespace Kratos.
