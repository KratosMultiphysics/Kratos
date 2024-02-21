//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/one-phase_flow/line_load_2D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition() : GeneralUPlDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPlDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineLoad2DDiffOrderCondition::LineLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPlDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineLoad2DDiffOrderCondition::~LineLoad2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineLoad2DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad2DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineLoad2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    Vector LineLoad = ZeroVector(3);
    rVariables.ConditionVector.resize(2,false);
    noalias(rVariables.ConditionVector) = ZeroVector(2);

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        LineLoad = rGeom[i].FastGetSolutionStepValue(FACE_LOAD);

        rVariables.ConditionVector[0] += rVariables.Nu[i]*LineLoad[0];
        rVariables.ConditionVector[1] += rVariables.Nu[i]*LineLoad[1];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoad2DDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rVariables.IntegrationCoefficient = ds * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoad2DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    SizeType Index;

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        Index = i * 2;

        rRightHandSideVector[Index]   += rVariables.Nu[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index+1] += rVariables.Nu[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
