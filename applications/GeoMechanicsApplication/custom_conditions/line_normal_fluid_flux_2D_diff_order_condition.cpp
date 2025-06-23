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
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

// Default Constructor
LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition()
    : LineLoad2DDiffOrderCondition()
{
}

// Constructor 1
LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition(IndexType NewId,
                                                                                 GeometryType::Pointer pGeometry)
    : LineLoad2DDiffOrderCondition(NewId, pGeometry)
{
}

// Constructor 2
LineNormalFluidFlux2DDiffOrderCondition::LineNormalFluidFlux2DDiffOrderCondition(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LineLoad2DDiffOrderCondition(NewId, pGeometry, pProperties)
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
    Vector nodal_normal_fluid_flux_vector(mpPressureGeometry->PointsNumber());
    VariablesUtilities::GetNodalValues(*mpPressureGeometry, NORMAL_FLUID_FLUX,
                                       nodal_normal_fluid_flux_vector.begin());
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

std::string LineNormalFluidFlux2DDiffOrderCondition::Info() const
{
    return "LineNormalFluidFlux2DDiffOrderCondition";
}

} // Namespace Kratos.
