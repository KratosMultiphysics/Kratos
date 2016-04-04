//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/surface_normal_fluid_flux_3D_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceNormalFluidFlux3DCondition::SurfaceNormalFluidFlux3DCondition() : GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalFluidFlux3DCondition::SurfaceNormalFluidFlux3DCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalFluidFlux3DCondition::SurfaceNormalFluidFlux3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalFluidFlux3DCondition::~SurfaceNormalFluidFlux3DCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalFluidFlux3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalFluidFlux3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    rVariables.ConditionVector.resize(1,false);
    rVariables.ConditionVector[0] = 0.0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rVariables.ConditionVector[0] += rVariables.Np[i]*rGeom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
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

void SurfaceNormalFluidFlux3DCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * (3 + 1) + 3;

        rRightHandSideVector[Global_i] -= rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
