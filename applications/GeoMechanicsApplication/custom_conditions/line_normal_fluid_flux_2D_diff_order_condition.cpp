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
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
LineNormalFluidFlux2DDiffOrderCondition::
    LineNormalFluidFlux2DDiffOrderCondition() : LineLoad2DDiffOrderCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
LineNormalFluidFlux2DDiffOrderCondition::
    LineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                            GeometryType::Pointer pGeometry) :
                                            LineLoad2DDiffOrderCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
LineNormalFluidFlux2DDiffOrderCondition::
    LineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                            GeometryType::Pointer pGeometry,
                                            PropertiesType::Pointer pProperties) :
                                            LineLoad2DDiffOrderCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer LineNormalFluidFlux2DDiffOrderCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineNormalFluidFlux2DDiffOrderCondition(NewId,
                                                                          GetGeometry().Create(ThisNodes),
                                                                          pProperties) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LineNormalFluidFlux2DDiffOrderCondition::
    CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    rVariables.ConditionVector.resize(1,false);
    rVariables.ConditionVector[0] = 0.0;

    for ( SizeType i = 0; i < NumPNodes; ++i ) {
        rVariables.ConditionVector[0] += rVariables.Np[i]*GetGeometry()[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void LineNormalFluidFlux2DDiffOrderCondition::
    CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                  ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes*2+i] -=  rVariables.Np[i]
                                              * rVariables.ConditionVector[0]
                                              * rVariables.IntegrationCoefficient;
    }
}

} // Namespace Kratos.
