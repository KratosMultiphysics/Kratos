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
#include "custom_utilities/variables_utilities.hpp"

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
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer SurfaceNormalFluidFlux3DDiffOrderCondition::Create(IndexType             NewId,
                                                                      GeometryType::Pointer pGeom,
                                                                      PropertiesType::Pointer pProperties) const
{
    return make_intrusive<SurfaceNormalFluidFlux3DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                          unsigned int PointNumber)
{
    KRATOS_TRY
    Vector nodal_normal_fluid_flux_vector(mpPressureGeometry->PointsNumber());
    VariablesUtilities::GetNodalValues(*mpPressureGeometry, NORMAL_FLUID_FLUX,
                                       nodal_normal_fluid_flux_vector.begin());
    rVariables.ConditionVector.resize(1, false);
    rVariables.ConditionVector[0] = MathUtils<>::Dot(rVariables.Np, nodal_normal_fluid_flux_vector);
    KRATOS_CATCH("")
}

void SurfaceNormalFluidFlux3DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                               ConditionVariables& rVariables)
{
    const SizeType num_u_nodes = GetGeometry().PointsNumber();
    for (SizeType node = 0; node < mpPressureGeometry->PointsNumber(); ++node) {
        rRightHandSideVector[num_u_nodes * 3 + node] -=
            rVariables.Np[node] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

std::string SurfaceNormalFluidFlux3DDiffOrderCondition::Info() const
{
    return "SurfaceNormalFluidFlux3DDiffOrderCondition";
}

} // Namespace Kratos.
