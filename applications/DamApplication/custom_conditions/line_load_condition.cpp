//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:               $
//

// Project includes
#include "custom_conditions/line_load_condition.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
LineLoadCondition::LineLoadCondition() : GeneralCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineLoadCondition::LineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineLoadCondition::LineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineLoadCondition::~LineLoadCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineLoadCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Vector LineLoad = ZeroVector(3);
    rVariables.ConditionVector = ZeroVector(2);

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        LineLoad = GetGeometry()[i].FastGetSolutionStepValue(LINE_LOAD);

        rVariables.ConditionVector[0] += rVariables.N[i]*LineLoad[0];
        rVariables.ConditionVector[1] += rVariables.N[i]*LineLoad[1];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoadCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rVariables.IntegrationCoefficient = GetProperties()[THICKNESS] * ds * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoadCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * 2;

        rRightHandSideVector[Global_i]   += rVariables.N[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+1] += rVariables.N[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }

}

} // Namespace Kratos.
