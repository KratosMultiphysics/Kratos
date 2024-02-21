//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/one-phase_flow/line_normal_liquid_flux_2D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
LineNormalLiquidFlux2DDiffOrderCondition::LineNormalLiquidFlux2DDiffOrderCondition() : GeneralUPlDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineNormalLiquidFlux2DDiffOrderCondition::LineNormalLiquidFlux2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPlDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineNormalLiquidFlux2DDiffOrderCondition::LineNormalLiquidFlux2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPlDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineNormalLiquidFlux2DDiffOrderCondition::~LineNormalLiquidFlux2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineNormalLiquidFlux2DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineNormalLiquidFlux2DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineNormalLiquidFlux2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    rVariables.ConditionVector.resize(1,false);
    rVariables.ConditionVector[0] = 0.0;

    for ( SizeType i = 0; i < NumPNodes; i++ )
    {
        rVariables.ConditionVector[0] += rVariables.Np[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_LIQUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLiquidFlux2DDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rVariables.IntegrationCoefficient = ds * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLiquidFlux2DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*2+i] -= rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
