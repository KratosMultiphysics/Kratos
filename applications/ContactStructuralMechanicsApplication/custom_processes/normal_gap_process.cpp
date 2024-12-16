// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "custom_processes/normal_gap_process.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>::Execute()
{
    KRATOS_TRY

    // We get the process info
    const ProcessInfo& r_process_info = mrMasterModelPart.GetRootModelPart().GetProcessInfo();

    // Iterate over the nodes
    auto& r_nodes_array_master = mrMasterModelPart.Nodes();
    auto& r_nodes_array_slave = mrSlaveModelPart.Nodes();

    // We set the auxiliary Coordinates
    const array_1d<double, 3> zero_array = ZeroVector(3);
    block_for_each(r_nodes_array_master, [&](Node& rNode) {
        if (mSearchOrientation) {
            rNode.SetValue(AUXILIAR_COORDINATES, rNode.Coordinates());
        } else {
            rNode.SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    });
    block_for_each(r_nodes_array_slave, [&](Node& rNode) {
        if (!mSearchOrientation) {
            rNode.SetValue(AUXILIAR_COORDINATES, rNode.Coordinates());
        } else {
            rNode.SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    });

    // Switch MASTER/SLAVE
    if (!mSearchOrientation) {
        SwitchFlagNodes(r_nodes_array_master);
        SwitchFlagNodes(r_nodes_array_slave);
    }

    // We set the mapper parameters
    Parameters mapping_parameters = Parameters(R"({"distance_threshold" : 1.0e24,"update_interface" : false, "remove_isolated_conditions" : true, "origin_variable_historical" : false, "destination_variable_historical" : false, "zero_tolerance_factor" : 1.0e0, "consider_tessellation" : false})" );
    if (r_process_info.Has(DISTANCE_THRESHOLD)) {
        mapping_parameters["distance_threshold"].SetDouble(r_process_info[DISTANCE_THRESHOLD]);
    }
    if (r_process_info.Has(ZERO_TOLERANCE_FACTOR)) {
        mapping_parameters["zero_tolerance_factor"].SetDouble(r_process_info[ZERO_TOLERANCE_FACTOR]);
    }
    const auto& r_properties = mrSlaveModelPart.Conditions().begin()->GetProperties();
    if (r_properties.Has(CONSIDER_TESSELLATION)) {
        mapping_parameters["consider_tessellation"].SetBool(r_properties[CONSIDER_TESSELLATION]);
    }
    MapperType mapper(mrMasterModelPart, mrSlaveModelPart, AUXILIAR_COORDINATES, mapping_parameters);
    mapper.Execute();

    // Switch again MASTER/SLAVE
    if (!mSearchOrientation) {
        SwitchFlagNodes(r_nodes_array_master);
        SwitchFlagNodes(r_nodes_array_slave);
    }

    // We compute now the normal gap and set the nodes under certain threshold as active
    ComputeNormalGap(r_nodes_array_master);
    ComputeNormalGap(r_nodes_array_slave);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>::ComputeNormalGap(NodesArrayType& rNodes)
{
    KRATOS_TRY

    struct auxiliary {array_1d<double, 3> normal, auxiliary_coordinates, components_gap; double gap = 0.0; };
    block_for_each(rNodes, auxiliary(), [this](Node& rNode, auxiliary& aux) {
        if (rNode.Is(SLAVE) == this->mSearchOrientation) {
            // We compute the gap
            noalias(aux.normal) = rNode.FastGetSolutionStepValue(NORMAL);
            noalias(aux.auxiliary_coordinates) = rNode.GetValue(AUXILIAR_COORDINATES);
            noalias(aux.components_gap) = ( rNode.Coordinates() - aux.auxiliary_coordinates);
            aux.gap = inner_prod(aux.components_gap, - aux.normal);

            // We activate if the node is close enough
            if (norm_2(aux.auxiliary_coordinates) > ZeroTolerance)
                rNode.SetValue(NORMAL_GAP, aux.gap);
        } else {
            rNode.SetValue(NORMAL_GAP, 0.0);
        }
    });

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class NormalGapProcess<2, 2>;
template class NormalGapProcess<3, 3>;
template class NormalGapProcess<3, 4>;
template class NormalGapProcess<3, 3, 4>;
template class NormalGapProcess<3, 4, 3>;

}  // namespace Kratos.
