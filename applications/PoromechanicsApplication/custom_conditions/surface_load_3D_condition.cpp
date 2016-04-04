//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/surface_load_3D_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceLoad3DCondition::SurfaceLoad3DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceLoad3DCondition::~SurfaceLoad3DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceLoad3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceLoad3DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    Vector SurfaceLoad = ZeroVector(3);
    rVariables.ConditionVector.resize(3,false);
    noalias(rVariables.ConditionVector) = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        SurfaceLoad = rGeom[i].FastGetSolutionStepValue(SURFACE_LOAD);

        rVariables.ConditionVector[0] += rVariables.Np[i]*SurfaceLoad[0];
        rVariables.ConditionVector[1] += rVariables.Np[i]*SurfaceLoad[1];
        rVariables.ConditionVector[2] += rVariables.Np[i]*SurfaceLoad[2];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceLoad3DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
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

void SurfaceLoad3DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
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
