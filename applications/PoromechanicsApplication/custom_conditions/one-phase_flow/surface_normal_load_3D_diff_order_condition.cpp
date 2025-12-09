//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/one-phase_flow/surface_normal_load_3D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition() : GeneralUPlDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : GeneralUPlDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalLoad3DDiffOrderCondition::SurfaceNormalLoad3DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : GeneralUPlDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalLoad3DDiffOrderCondition::~SurfaceNormalLoad3DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalLoad3DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalLoad3DDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalLoad3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    double NormalVector[3];

    NormalVector[0] = rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](2,1) -
                      rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](1,1);

    NormalVector[1] = rVariables.JContainer[PointNumber](2,0) * rVariables.JContainer[PointNumber](0,1) -
                      rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](2,1);

    NormalVector[2] = rVariables.JContainer[PointNumber](0,0) * rVariables.JContainer[PointNumber](1,1) -
                      rVariables.JContainer[PointNumber](1,0) * rVariables.JContainer[PointNumber](0,1);
    
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    double NormalStress = 0.0;
    rVariables.ConditionVector.resize(3,false);

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        NormalStress += rVariables.Nu[i]*rGeom[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    }

    rVariables.ConditionVector[0] = NormalStress * NormalVector[0];
    rVariables.ConditionVector[1] = NormalStress * NormalVector[1];
    rVariables.ConditionVector[2] = NormalStress * NormalVector[2];

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoad3DDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    rVariables.IntegrationCoefficient = weight;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SurfaceNormalLoad3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    SizeType Index;

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        Index = i * 3;
        
        rRightHandSideVector[Index]   += rVariables.Nu[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index+1] += rVariables.Nu[i] * rVariables.ConditionVector[1] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Index+2] += rVariables.Nu[i] * rVariables.ConditionVector[2] * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
