// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
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
    const auto it_node_begin_master = r_nodes_array_master.begin();
    auto& r_nodes_array_slave = mrSlaveModelPart.Nodes();
    const auto it_node_begin_slave = r_nodes_array_slave.begin();

    // We set the auxiliar Coordinates
    const array_1d<double, 3> zero_array = ZeroVector(3);
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array_master.size()); ++i) {
        auto it_node = it_node_begin_master + i;

        if (mSearchOrientation) {
            it_node->SetValue(AUXILIAR_COORDINATES, it_node->Coordinates());
        } else {
            it_node->SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array_slave.size()); ++i) {
        auto it_node = it_node_begin_slave + i;

        if (!mSearchOrientation) {
            it_node->SetValue(AUXILIAR_COORDINATES, it_node->Coordinates());
        } else {
            it_node->SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    }

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

    array_1d<double, 3> normal, auxiliar_coordinates, components_gap;
    double gap = 0.0;

    // The first iterator
    const auto it_node_begin = rNodes.begin();

    #pragma omp parallel for firstprivate(gap, normal, auxiliar_coordinates, components_gap)
    for(int i = 0; i < static_cast<int>(rNodes.size()); ++i) {
        auto it_node = it_node_begin + i;

        if (it_node->Is(SLAVE) == mSearchOrientation) {
            // We compute the gap
            noalias(normal) = it_node->FastGetSolutionStepValue(NORMAL);
            noalias(auxiliar_coordinates) = it_node->GetValue(AUXILIAR_COORDINATES);
            noalias(components_gap) = ( it_node->Coordinates() - auxiliar_coordinates);
            gap = inner_prod(components_gap, - normal);

            // We activate if the node is close enough
            if (norm_2(auxiliar_coordinates) > ZeroTolerance)
                it_node->SetValue(NORMAL_GAP, gap);
        } else {
            it_node->SetValue(NORMAL_GAP, 0.0);
        }
    }

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
