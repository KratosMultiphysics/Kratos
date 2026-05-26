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
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.h"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition()
    : SurfaceLoad3DDiffOrderCondition()
{
}

SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(IndexType NewId,
                                                                                       GeometryType::Pointer pGeometry)
    : SurfaceLoad3DDiffOrderCondition(NewId, std::move(pGeometry))
{
}

SurfaceNormalFluidFlux3DDiffOrderCondition::SurfaceNormalFluidFlux3DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SurfaceLoad3DDiffOrderCondition(NewId, std::move(pGeometry), std::move(pProperties))
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
    const auto nodal_normal_fluid_flux_vector =
        VariablesUtilities::GetNodalValues(*mpPressureGeometry, NORMAL_FLUID_FLUX);
    rVariables.ConditionVector =
        ScalarVector{1, std::inner_product(rVariables.Np.cbegin(), rVariables.Np.cend(),
                                           nodal_normal_fluid_flux_vector.cbegin(), 0.0)};
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

void SurfaceNormalFluidFlux3DDiffOrderCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SurfaceLoad3DDiffOrderCondition)
}

void SurfaceNormalFluidFlux3DDiffOrderCondition::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SurfaceLoad3DDiffOrderCondition)}

std::string SurfaceNormalFluidFlux3DDiffOrderCondition::Info() const
{
    return "SurfaceNormalFluidFlux3DDiffOrderCondition";
}

} // Namespace Kratos.
