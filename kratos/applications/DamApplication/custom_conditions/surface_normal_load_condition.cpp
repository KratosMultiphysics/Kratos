//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:        $
//

// Project includes
#include "custom_conditions/surface_normal_load_condition.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
SurfaceNormalLoadCondition::SurfaceNormalLoadCondition() : GeneralCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalLoadCondition::SurfaceNormalLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalLoadCondition::SurfaceNormalLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalLoadCondition::~SurfaceNormalLoadCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalLoadCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    double NormalVector[3];

    NormalVector[0] = rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](2,1) -
                      rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](1,1);

    NormalVector[1] = rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](0,1) -
                      rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](2,1);

    NormalVector[2] = rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](1,1) -
                      rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](0,1);

    const unsigned int number_of_nodes = GetGeometry().size();
    double NormalStress = 0;
    rVariables.ConditionVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        NormalStress += rVariables.N[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    }

    rVariables.ConditionVector[0] = NormalStress * NormalVector[0];
    rVariables.ConditionVector[1] = NormalStress * NormalVector[1];
    rVariables.ConditionVector[2] = NormalStress * NormalVector[2];


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoadCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    rVariables.IntegrationCoefficient = weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoadCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
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
