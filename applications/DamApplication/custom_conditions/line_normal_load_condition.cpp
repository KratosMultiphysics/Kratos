//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:        $
//

// Project includes
#include "custom_conditions/line_normal_load_condition.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
LineNormalLoadCondition::LineNormalLoadCondition() : GeneralCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineNormalLoadCondition::LineNormalLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineNormalLoadCondition::LineNormalLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineNormalLoadCondition::~LineNormalLoadCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineNormalLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineNormalLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineNormalLoadCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    double NormalStress = 0;
    double TangentialStress = 0;
    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);
    rVariables.ConditionVector.resize(2,false);
    noalias(rVariables.ConditionVector) = ZeroVector(2);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        NormalStress     += rVariables.N[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
        TangentialStress += rVariables.N[i]*GetGeometry()[i].FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
    }

    rVariables.ConditionVector[0] = TangentialStress * dx_dxi - NormalStress     * dy_dxi;
    rVariables.ConditionVector[1] = NormalStress     * dx_dxi + TangentialStress * dy_dxi;


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLoadCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    rVariables.IntegrationCoefficient = GetProperties()[THICKNESS] * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLoadCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
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
