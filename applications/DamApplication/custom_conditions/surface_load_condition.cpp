//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:        $
//

// Project includes
#include "custom_conditions/surface_load_condition.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
SurfaceLoadCondition::SurfaceLoadCondition() : GeneralCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceLoadCondition::SurfaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceLoadCondition::SurfaceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceLoadCondition::~SurfaceLoadCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceLoadCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Vector SurfaceLoad = ZeroVector(3);
    rVariables.ConditionVector.resize(3);
    rVariables.ConditionVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        SurfaceLoad = GetGeometry()[i].FastGetSolutionStepValue(SURFACE_LOAD);

        rVariables.ConditionVector[0] += rVariables.N[i]*SurfaceLoad[0];
        rVariables.ConditionVector[1] += rVariables.N[i]*SurfaceLoad[1];
        rVariables.ConditionVector[2] += rVariables.N[i]*SurfaceLoad[2];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceLoadCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    double NormalVector[3];

    NormalVector[0] = rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](2,1) -
                      rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](1,1);

    NormalVector[1] = rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](0,1) -
                      rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](2,1);

    NormalVector[2] = rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](1,1) -
                      rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](0,1);

    double dA = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);

    rVariables.IntegrationCoefficient = dA * weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceLoadCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * 3;

        rRightHandSideVector[Global_i]   += rVariables.N[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+1] += rVariables.N[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+2] += rVariables.N[i] * rVariables.ConditionVector[2] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
