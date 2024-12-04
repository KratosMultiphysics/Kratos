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
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition() : GeneralUPwDiffOrderCondition()
{
}

// Constructor 1
SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeneralUPwDiffOrderCondition(NewId, pGeometry)
{
}

// Constructor 2
SurfaceLoad3DDiffOrderCondition::SurfaceLoad3DDiffOrderCondition(IndexType             NewId,
                                                                 GeometryType::Pointer pGeometry,
                                                                 PropertiesType::Pointer pProperties)
    : GeneralUPwDiffOrderCondition(NewId, pGeometry, pProperties)
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

    const GeometryType& rGeom       = GetGeometry();
    const SizeType      NumUNodes   = rGeom.PointsNumber();
    Vector              SurfaceLoad = ZeroVector(3);
    rVariables.ConditionVector.resize(3, false);
    noalias(rVariables.ConditionVector) = ZeroVector(3);

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SurfaceLoad = rGeom[i].FastGetSolutionStepValue(SURFACE_LOAD);

        rVariables.ConditionVector[0] += rVariables.Nu[i] * SurfaceLoad[0];
        rVariables.ConditionVector[1] += rVariables.Nu[i] * SurfaceLoad[1];
        rVariables.ConditionVector[2] += rVariables.Nu[i] * SurfaceLoad[2];
    }

    KRATOS_CATCH("")
}

double SurfaceLoad3DDiffOrderCondition::CalculateIntegrationCoefficient(
    const IndexType                                 PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    double NormalVector[3];

    NormalVector[0] = JContainer[PointNumber](1, 0) * JContainer[PointNumber](2, 1) -
                      JContainer[PointNumber](2, 0) * JContainer[PointNumber](1, 1);

    NormalVector[1] = JContainer[PointNumber](2, 0) * JContainer[PointNumber](0, 1) -
                      JContainer[PointNumber](0, 0) * JContainer[PointNumber](2, 1);

    NormalVector[2] = JContainer[PointNumber](0, 0) * JContainer[PointNumber](1, 1) -
                      JContainer[PointNumber](1, 0) * JContainer[PointNumber](0, 1);

    double dA = sqrt(NormalVector[0] * NormalVector[0] + NormalVector[1] * NormalVector[1] +
                     NormalVector[2] * NormalVector[2]);

    return dA * IntegrationPoints[PointNumber].Weight();

    KRATOS_CATCH("")
}

void SurfaceLoad3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
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

std::string SurfaceLoad3DDiffOrderCondition::Info() const
{
    return "SurfaceLoad3DDiffOrderCondition";
}

} // Namespace Kratos.
