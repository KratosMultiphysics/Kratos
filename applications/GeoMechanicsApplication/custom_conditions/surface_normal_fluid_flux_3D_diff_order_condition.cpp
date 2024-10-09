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
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition()
    : SurfaceLoad3DDiffOrderCondition()
{
}

// Constructor 1
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(IndexType NewId,
                                                                                       GeometryType::Pointer pGeometry)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry)
{
}

// Constructor 2
SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SurfaceLoad3DDiffOrderCondition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer SurfaceNormalFluidFlux3DDiffOrderCondition::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalFluidFlux3DDiffOrderCondition(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                          unsigned int PointNumber)
{
    KRATOS_TRY

    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    rVariables.ConditionVector.resize(1, false);
    rVariables.ConditionVector[0] = 0.0;

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rVariables.ConditionVector[0] +=
            rVariables.Np[i] * GetGeometry()[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }

    KRATOS_CATCH("")
}

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector,
                                                                               ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes * 3 + i] -=
            rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

std::string SurfaceNormalFluidFlux3DDiffOrderCondition::Info() const
{
    return "SurfaceNormalFluidFlux3DDiffOrderCondition";
}

} // Namespace Kratos.
