//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/one-phase_flow/surface_normal_liquid_flux_3D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceNormalLiquidFlux3DDiffOrderCondition::SurfaceNormalLiquidFlux3DDiffOrderCondition() : GeneralUPlDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalLiquidFlux3DDiffOrderCondition::SurfaceNormalLiquidFlux3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPlDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalLiquidFlux3DDiffOrderCondition::SurfaceNormalLiquidFlux3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPlDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalLiquidFlux3DDiffOrderCondition::~SurfaceNormalLiquidFlux3DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalLiquidFlux3DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalLiquidFlux3DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalLiquidFlux3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    rVariables.ConditionVector.resize(1,false);
    rVariables.ConditionVector[0] = 0.0;

    for ( SizeType i = 0; i < NumPNodes; i++ )
    {
        rVariables.ConditionVector[0] += rVariables.Np[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_LIQUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLiquidFlux3DDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
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

void SurfaceNormalLiquidFlux3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rRightHandSideVector[NumUNodes*3+i] -= rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
