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
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition() : GeneralUPwDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

// Constructor 1
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeneralUPwDiffOrderCondition(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------

// Constructor 2
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType               NewId,
                                                           GeometryType::Pointer   pGeometry,
                                                           PropertiesType::Pointer pProperties)
    : GeneralUPwDiffOrderCondition(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineLoad2DDiffOrderCondition::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad2DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineLoad2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const GeometryType& rGeom     = GetGeometry();
    const SizeType      NumUNodes = rGeom.PointsNumber();
    Vector              LineLoad  = ZeroVector(3);
    rVariables.ConditionVector.resize(2, false);
    noalias(rVariables.ConditionVector) = ZeroVector(2);

    for (SizeType i = 0; i < NumUNodes; ++i) {
        LineLoad = rGeom[i].FastGetSolutionStepValue(LINE_LOAD);

        rVariables.ConditionVector[0] += rVariables.Nu[i] * LineLoad[0];
        rVariables.ConditionVector[1] += rVariables.Nu[i] * LineLoad[1];
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
double LineLoad2DDiffOrderCondition::CalculateIntegrationCoefficient(
    const IndexType                                 PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    const double dx_dxi = JContainer[PointNumber](0, 0);
    const double dy_dxi = JContainer[PointNumber](1, 0);

    const double ds = sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

    return ds * IntegrationPoints[PointNumber].Weight();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void LineLoad2DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                                                 ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SizeType Index = i * 2;

        rRightHandSideVector[Index] +=
            rVariables.Nu[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index + 1] +=
            rVariables.Nu[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
