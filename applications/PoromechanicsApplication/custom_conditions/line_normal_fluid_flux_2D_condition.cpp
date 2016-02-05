//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/line_normal_fluid_flux_2D_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
LineNormalFluidFlux2DCondition::LineNormalFluidFlux2DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineNormalFluidFlux2DCondition::LineNormalFluidFlux2DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineNormalFluidFlux2DCondition::LineNormalFluidFlux2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
LineNormalFluidFlux2DCondition::~LineNormalFluidFlux2DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineNormalFluidFlux2DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineNormalFluidFlux2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineNormalFluidFlux2DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    rVariables.ConditionVector.resize(1,false);
    noalias(rVariables.ConditionVector) = ZeroVector(1);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rVariables.ConditionVector[0] += rVariables.Np[i]*rGeom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalFluidFlux2DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double dx_dxi = rVariables.JContainer[PointNumber](0,0), dy_dxi = rVariables.JContainer[PointNumber](1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rVariables.IntegrationCoefficient = GetProperties()[THICKNESS] * ds * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalFluidFlux2DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * (2 + 1) + 2;

        rRightHandSideVector[Global_i] -= rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
