//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/line_normal_load_2D_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
LineNormalLoad2DCondition::LineNormalLoad2DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineNormalLoad2DCondition::LineNormalLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineNormalLoad2DCondition::LineNormalLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineNormalLoad2DCondition::~LineNormalLoad2DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineNormalLoad2DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineNormalLoad2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineNormalLoad2DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    double NormalStress = 0.0;
    double TangentialStress = 0.0;
    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);
    rVariables.ConditionVector.resize(2,false);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        NormalStress     += rVariables.Np[i]*rGeom[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
        TangentialStress += rVariables.Np[i]*rGeom[i].FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
    }

    rVariables.ConditionVector[0] = TangentialStress * dx_dxi - NormalStress     * dy_dxi;
    rVariables.ConditionVector[1] = NormalStress     * dx_dxi + TangentialStress * dy_dxi;


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLoad2DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    rVariables.IntegrationCoefficient = GetProperties()[THICKNESS] * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalLoad2DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
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
