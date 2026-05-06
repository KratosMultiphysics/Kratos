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
#include "custom_conditions/line_load_2D_diff_order_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"

namespace Kratos
{

LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition() : GeneralUPwDiffOrderCondition() {}

LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeneralUPwDiffOrderCondition(NewId, std::move(pGeometry))
{
}

LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType               NewId,
                                                           GeometryType::Pointer   pGeometry,
                                                           PropertiesType::Pointer pProperties)
    : GeneralUPwDiffOrderCondition(NewId, std::move(pGeometry), std::move(pProperties))
{
}

Condition::Pointer LineLoad2DDiffOrderCondition::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LineLoad2DDiffOrderCondition::Create(IndexType               NewId,
                                                        GeometryType::Pointer   pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return make_intrusive<LineLoad2DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void LineLoad2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const auto line_load_vectors = VariablesUtilities::GetNodalValues(GetGeometry(), LINE_LOAD);

    rVariables.ConditionVector.resize(2, false);
    noalias(rVariables.ConditionVector) = ZeroVector(2);
    for (SizeType i = 0; i < line_load_vectors.size(); ++i) {
        rVariables.ConditionVector[0] += rVariables.Nu[i] * line_load_vectors[i][0];
        rVariables.ConditionVector[1] += rVariables.Nu[i] * line_load_vectors[i][1];
    }

    KRATOS_CATCH("")
}

double LineLoad2DDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    return ConditionUtilities::CalculateIntegrationCoefficient(
        JContainer[PointNumber], IntegrationPoints[PointNumber].Weight());
}

void LineLoad2DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                 ConditionVariables& rVariables)
{
    for (SizeType node = 0; node < this->GetGeometry().PointsNumber(); ++node) {
        rRightHandSideVector[2 * node] +=
            rVariables.Nu[node] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[2 * node + 1] +=
            rVariables.Nu[node] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }
}

void LineLoad2DDiffOrderCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeneralUPwDiffOrderCondition)
}

void LineLoad2DDiffOrderCondition::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeneralUPwDiffOrderCondition)}

std::string LineLoad2DDiffOrderCondition::Info() const
{
    return "LineLoad2DDiffOrderCondition";
}

} // Namespace Kratos.
