//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_conditions/line_load_2D_condition.hpp"

namespace Kratos
{

// Default Constructor
LineLoad2DCondition::LineLoad2DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineLoad2DCondition::LineLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineLoad2DCondition::LineLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineLoad2DCondition::~LineLoad2DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineLoad2DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineLoad2DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    Vector LineLoad = ZeroVector(3);
    rVariables.ConditionVector.resize(2,false);
    noalias(rVariables.ConditionVector) = ZeroVector(2);

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        LineLoad = rGeom[i].FastGetSolutionStepValue(LINE_LOAD);

        rVariables.ConditionVector[0] += rVariables.Np[i]*LineLoad[0];
        rVariables.ConditionVector[1] += rVariables.Np[i]*LineLoad[1];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoad2DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rVariables.IntegrationCoefficient = GetProperties()[THICKNESS] * ds * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineLoad2DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * (2 + 1);

        rRightHandSideVector[Global_i]   += rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+1] += rVariables.Np[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
    }

}

} // Namespace Kratos.
