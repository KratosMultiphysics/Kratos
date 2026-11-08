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
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.h"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition()
    : LineLoad2DDiffOrderCondition()
{
}

LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                                                                 GeometryType::Pointer pGeometry)
    : LineLoad2DDiffOrderCondition(NewId, std::move(pGeometry))
{
}

LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LineLoad2DDiffOrderCondition(NewId, std::move(pGeometry), std::move(pProperties))
{
}

Condition::Pointer LineNormalFluidFlux2DDiffOrderCondition::Create(IndexType             NewId,
                                                                   NodesArrayType const& ThisNodes,
                                                                   PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LineNormalFluidFlux2DDiffOrderCondition::Create(IndexType             NewId,
                                                                   GeometryType::Pointer pGeom,
                                                                   PropertiesType::Pointer pProperties) const
{
    return make_intrusive<LineNormalFluidFlux2DDiffOrderCondition>(NewId, pGeom, pProperties);
}

void LineNormalFluidFlux2DDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables,
                                                                       unsigned int PointNumber)
{
    KRATOS_TRY
    const auto nodal_normal_fluid_flux_vector =
        VariablesUtilities::GetNodalValues(*mpPressureGeometry, NORMAL_FLUID_FLUX);
    rVariables.ConditionVector =
        ScalarVector(1, std::inner_product(rVariables.Np.cbegin(), rVariables.Np.cend(),
                                           nodal_normal_fluid_flux_vector.cbegin(), 0.0));
    KRATOS_CATCH("")
}

void LineNormalFluidFlux2DDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                            ConditionVariables& rVariables)
{
    const SizeType NumUNodes = GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes * 2 + i] -=
            rVariables.Np[i] * rVariables.ConditionVector[0] * rVariables.IntegrationCoefficient;
    }
}

void LineNormalFluidFlux2DDiffOrderCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LineLoad2DDiffOrderCondition)
}

void LineNormalFluidFlux2DDiffOrderCondition::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LineLoad2DDiffOrderCondition)}

std::string LineNormalFluidFlux2DDiffOrderCondition::Info() const
{
    return "LineNormalFluidFlux2DDiffOrderCondition";
}

} // Namespace Kratos.
