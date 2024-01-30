//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah,
//                   Ruben Zorilla
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "wss_statistics_utilities.h"
#include "fluid_dynamics_biomedical_application_variables.h"

namespace Kratos {

void WssStatisticsUtilities::InitializeWSSVariables(ModelPart& rModelPart)
{
    // Initialize WSS variables
    const array_1d<double,3> aux_zero = ZeroVector(3);
    block_for_each(rModelPart.Nodes(), [&aux_zero](NodeType& rNode){
        rNode.SetValue(WSS, 0.0);
        rNode.SetValue(RRT, 0.0);
        rNode.SetValue(OSI, 0.0);
        rNode.SetValue(ECAP, 0.0);
        rNode.SetValue(TWSS, 0.0);
        rNode.SetValue(TAWSS, 0.0);
        rNode.SetValue(TEMPORAL_OSI, aux_zero);
        rNode.SetValue(WSS_NORMAL_STRESS, aux_zero);
        rNode.SetValue(WSS_TANGENTIAL_STRESS, aux_zero);
    });
}

void WssStatisticsUtilities::CalculateWSS(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rNormalVariable,
    const bool IsNormalHistorical)
{
    // Set the nodal normal getter
    std::function<const array_1d<double,3>(const NodeType&)> normal_getter;
    if (IsNormalHistorical) {
        normal_getter = [&rNormalVariable](const NodeType& rNode)->const array_1d<double,3>{return rNode.FastGetSolutionStepValue(rNormalVariable);};
    } else {
        normal_getter = [&rNormalVariable](const NodeType& rNode)->const array_1d<double,3>{return rNode.GetValue(rNormalVariable);};
    }

    // Declare the auxiliary TLS container
    struct WSSTLS
    {
        array_1d<double,3> normal_tls;
        array_1d<double,3> normal_load;
        array_1d<double,3> tangential_load;
    };

    // Loop the WSS model part nodes
    block_for_each(rModelPart.Nodes(), WSSTLS(), [&normal_getter](NodeType& rNode, WSSTLS& rAuxTLS){
        // Get TLS arrays
        auto& normal = rAuxTLS.normal_tls;
        auto& normal_load = rAuxTLS.normal_load;
        auto& tangential_load = rAuxTLS.tangential_load;

        // Normalize nodal normal
        noalias(normal) = normal_getter(rNode);
        const double normal_norm = norm_2(normal);
        if (normal_norm > 1.0e-12) {
            normal /= normal_norm;
        } else {
            KRATOS_WARNING("CalculateWSS") << rNode.Id() << " has normal norm equal to " << normal_norm << "." << std::endl;
        }

        // Calculate the distributed projections
        const array_1d<double,3>& r_reaction = rNode.FastGetSolutionStepValue(REACTION);
        const double projection = inner_prod(r_reaction, normal);
        noalias(normal_load) = projection * normal;
        noalias(tangential_load) = r_reaction - normal_load;
        normal_load /= normal_norm;
        tangential_load /= normal_norm;

        // Save computed magnitudes
        rNode.GetValue(WSS) = norm_2(tangential_load);
        rNode.GetValue(WSS_NORMAL_STRESS) = normal_load;
        rNode.GetValue(WSS_TANGENTIAL_STRESS) = tangential_load;
    });
}

void WssStatisticsUtilities::CalculateOSI(ModelPart &rModelPart)
{
    const double abs_tol = std::numeric_limits<double>::epsilon();
    KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(TIME)) << "'TIME' is not present in '" << rModelPart.FullName() << "' ProcessInfo container." << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(DELTA_TIME)) << "'DELTA_TIME' is not present in '" << rModelPart.FullName() << "' ProcessInfo container." << std::endl;
    const double time = rModelPart.GetProcessInfo()[TIME];
    const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
    block_for_each(rModelPart.Nodes(), [&time, &delta_time, &abs_tol](NodeType& rNode){

        // Accumulate the current step contribution to the Oscillatory Shear Index (OSI)
        // Note that this requires the tangential WSS to be already computed (see CalculateWSS)
        auto& r_temporal_osi = rNode.GetValue(TEMPORAL_OSI);
        const auto& r_wss_tang_stress = rNode.GetValue(WSS_TANGENTIAL_STRESS);
        r_temporal_osi += r_wss_tang_stress*delta_time;

        // Calculates the sum of the WSS magnitudes for all time steps
        double& r_twss = rNode.GetValue(TWSS);
        double& r_tawss = rNode.GetValue(TAWSS);
        r_twss += norm_2(r_wss_tang_stress)*delta_time;
        r_tawss = r_twss / time;

        // Calculate OSI and OSI-dependent magnitudes
        double& r_osi = rNode.GetValue(OSI);
        const double temp_osi_norm = norm_2(r_temporal_osi);
        r_osi = temp_osi_norm / r_twss > 1.0 ? 0.0 : 0.5 * (1.0 - temp_osi_norm/r_twss);
        if (r_tawss > abs_tol) {
            rNode.GetValue(ECAP) = r_osi / r_tawss;
            if (std::abs(r_osi - 0.5) > abs_tol) {
                rNode.GetValue(RRT) = 1.0 / ((1.0 - 2.0*r_osi) * r_tawss);
            }
        }
    });
}

}
