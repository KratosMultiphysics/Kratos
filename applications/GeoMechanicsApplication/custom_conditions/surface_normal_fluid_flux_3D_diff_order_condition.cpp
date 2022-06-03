// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//


// Project includes
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition() : GeneralUPwDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPwDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPwDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalFluidFlux3DDiffOrderCondition::~SurfaceNormalFluidFlux3DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalFluidFlux3DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalFluidFlux3DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    rVariables.ConditionVector.resize(1,false);
    rVariables.ConditionVector[0] = 0.0;

    for ( SizeType i = 0; i < NumPNodes; i++ )
    {
        rVariables.ConditionVector[0] += rVariables.Np[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
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

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*3+i] -= rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
