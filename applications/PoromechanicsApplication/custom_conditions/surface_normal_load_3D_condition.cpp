//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/surface_normal_load_3D_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
SurfaceNormalLoad3DCondition::SurfaceNormalLoad3DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalLoad3DCondition::SurfaceNormalLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalLoad3DCondition::SurfaceNormalLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalLoad3DCondition::~SurfaceNormalLoad3DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalLoad3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalLoad3DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
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
        NormalStress += rVariables.Np[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    }

    rVariables.ConditionVector[0] = NormalStress * NormalVector[0];
    rVariables.ConditionVector[1] = NormalStress * NormalVector[1];
    rVariables.ConditionVector[2] = NormalStress * NormalVector[2];


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoad3DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    rVariables.IntegrationCoefficient = weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoad3DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * (3 + 1);

        rRightHandSideVector[Global_i]   += rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+1] += rVariables.Np[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+2] += rVariables.Np[i] * rVariables.ConditionVector[2] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
