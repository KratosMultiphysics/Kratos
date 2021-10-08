//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_conditions/data_containers/k_omega/condition_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_vis_log_law_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaWallConditionData
{
const Variable<double>& OmegaVisLogBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

void OmegaVisLogBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = rCondition.GetGeometry();
    const int number_of_nodes = r_geometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

void OmegaVisLogBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mWallDistance3 = std::pow(this->GetGeometry().GetValue(DISTANCE), 3);

    KRATOS_CATCH("");
}

double OmegaVisLogBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    const double g_n = 12.0 * rParameters.mKinematicViscosity / (0.0075 * mWallDistance3);

    return (rParameters.mKinematicViscosity) * g_n;
}

} // namespace KOmegaWallConditionData

} // namespace Kratos