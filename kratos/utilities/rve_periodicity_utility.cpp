//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "utilities/rve_periodicity_utility.h"
#include "constraints/linear_master_slave_constraint.h"

namespace Kratos
{
void RVEPeriodicityUtility::AssignPeriodicity(
    ModelPart& rMasterModelPart,
    ModelPart& rSlaveModelPart,
    const Matrix& rStrainTensor,
    const Vector& rDirection,
    const double SearchTolerance
    )
{
    // Get DOMAIN_SIZE from ProcessInfo container
    const auto& r_process_info = rMasterModelPart.GetProcessInfo();
    KRATOS_ERROR_IF_NOT(r_process_info.Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE not found in " << rMasterModelPart.FullName() << " model part ProcessInfo." << std::endl;
    const std::size_t domain_size = r_process_info[DOMAIN_SIZE];

    // Call the assign periodicity methods specialization
    if (domain_size == 2) {
        AuxiliaryAssignPeriodicity<2>(rMasterModelPart, rSlaveModelPart, rStrainTensor, rDirection, SearchTolerance);
    } else if (domain_size == 3) {
        AuxiliaryAssignPeriodicity<3>(rMasterModelPart, rSlaveModelPart, rStrainTensor, rDirection, SearchTolerance);
    } else {
        KRATOS_ERROR << "Wrong DOMAIN_SIZE value " << domain_size << " in " << rMasterModelPart.FullName() << " model part ProcessInfo." << std::endl;
    }
}

template<std::size_t TDim>
void RVEPeriodicityUtility::AuxiliaryAssignPeriodicity(
    ModelPart& rMasterModelPart,
    ModelPart& rSlaveModelPart,
    const Matrix& rStrainTensor,
    const Vector& rDirection,
    const double SearchTolerance
    )
{
    KRATOS_ERROR_IF(rMasterModelPart.NumberOfConditions() == 0) << "The master is expected to have conditions and it is empty" << std::endl;

    const Vector translation = prod(rStrainTensor, rDirection);
    array_1d<double,3> aux_direction = ZeroVector(3);
    for (std::size_t i = 0; i < rDirection.size(); ++i) {
        aux_direction(i) = rDirection[i];
    }

    BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(rMasterModelPart);
    bin_based_point_locator.UpdateSearchDatabase();

    int max_search_results = 100;

    // Construct auxiliary data structure to contain the master slave relation.
    // Slave nodes must appear once, however a non-circular dependency is allowed between the masters
    for (IndexType i = 0; i < rSlaveModelPart.Nodes().size(); ++i) {
        // Search in which condition it falls
        auto it_node = rSlaveModelPart.NodesBegin() + i;

        Condition::Pointer p_host_cond;
        Vector N;
        array_1d<double, 3> transformed_slave_coordinates = it_node->GetInitialPosition() - aux_direction;

        // Finding the host element for this node
        const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(transformed_slave_coordinates, N, p_host_cond, max_search_results, SearchTolerance);
        if (is_found) {
            const auto& r_geometry = p_host_cond->GetGeometry();

            DataTupletype aux_data;

            auto &T = std::get<2>(aux_data);
            T = translation;

            auto& r_master_ids = std::get<0>(aux_data);
            auto& r_weights = std::get<1>(aux_data);
            for (IndexType j = 0; j < r_geometry.size(); ++j) {
                r_master_ids.push_back(r_geometry[j].Id());
                r_weights.push_back(N[j]);
            }

            if (mAuxPairings.find(it_node->Id()) == mAuxPairings.end()) { // This slave is not already present
                mAuxPairings[it_node->Id()] = aux_data;
            } else {
                KRATOS_INFO("RVEPeriodicityUtility") << "Slave model part = " << rSlaveModelPart << std::endl;
                KRATOS_INFO("RVEPeriodicityUtility") << "Master model part = " << rMasterModelPart << std::endl;
                KRATOS_ERROR << "Attempting to add twice the slave node with Id " << it_node->Id() << std::endl;
            }
        } else {

            double node_distance = 1e50;
            IndexType closest_node_id = 0;

            DataTupletype aux_data;

            auto &T = std::get<2>(aux_data);
            T = translation;

            //brute force search amongst all nodes to find the closest
            for(auto& rMasterNode : rMasterModelPart.Nodes()) {
                array_1d<double,3> distance_vector = rMasterNode.GetInitialPosition() - transformed_slave_coordinates;
                double d = norm_2(distance_vector);
                if(d < node_distance) {
                    node_distance = d;
                    closest_node_id = rMasterNode.Id();
                }
            }

            auto& r_master_ids = std::get<0>(aux_data);
            auto& r_weights = std::get<1>(aux_data);
            r_master_ids.push_back(closest_node_id);
            r_weights.push_back(1.0);

            if (mAuxPairings.find(it_node->Id()) == mAuxPairings.end()) { // This slave is not already present
                mAuxPairings[it_node->Id()] = aux_data;
            } else {
                KRATOS_INFO("RVEPeriodicityUtility") << "Slave model part = " << rSlaveModelPart << std::endl;
                KRATOS_INFO("RVEPeriodicityUtility") << "Master model part = " << rMasterModelPart << std::endl;
                KRATOS_ERROR << "Attempting to add twice the slave node with Id " << it_node->Id() << std::endl;
            }

            std::cout << "brute force search assigned by to node " << it_node->Id() << " the node " << closest_node_id << " by fallback closest node search.  Distance is " <<  node_distance << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RVEPeriodicityUtility::AppendIdsAndWeights(
    std::map<IndexType, DataTupletype>& rAux,
    const IndexType MasterId,
    const double MasterWeight,
    std::vector<IndexType>& rFinalMastersIds,
    std::vector<double>& rFinalMastersWeights,
    Vector& rFinalT)
{
    if (std::abs(MasterWeight) > 1e-12) { // Discard nodes with negligible weight (note that weights sum to 1)
        if (rAux.find(MasterId) == rAux.end()) { // Master is NOT also a slave
            rFinalMastersIds.push_back(MasterId);
            rFinalMastersWeights.push_back(MasterWeight);
        } else { // Master also happens to be a slave
            const auto& r_other_data = rAux[MasterId];
            const auto& r_other_master_ids = std::get<0>(r_other_data);
            const auto& r_other_master_weights = std::get<1>(r_other_data);
            const auto& r_other_T = std::get<2>(r_other_data);
            for (IndexType j = 0; j < r_other_master_ids.size(); ++j) {
                AppendIdsAndWeights(rAux, r_other_master_ids[j], MasterWeight * r_other_master_weights[j], rFinalMastersIds, rFinalMastersWeights, rFinalT);
            }

            rFinalT += MasterWeight * r_other_T;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

MasterSlaveConstraint::Pointer RVEPeriodicityUtility::GenerateConstraint(
    IndexType& rConstraintId,
    const DoubleVariableType& rVar,
    NodeType::Pointer pSlaveNode,
    const std::vector<IndexType>& rMasterIds,
    const Matrix& rRelationMatrix,
    const Vector& rTranslationVector)
{
    DofPointerVectorType slave_dofs, master_dofs;
    slave_dofs.reserve(1);
    master_dofs.reserve(rMasterIds.size());

    slave_dofs.push_back(pSlaveNode->pGetDof(rVar));
    for (IndexType i = 0; i < rMasterIds.size(); ++i)
        master_dofs.push_back(mrModelPart.pGetNode(rMasterIds[i])->pGetDof(rVar));

    auto pconstraint = Kratos::make_shared<LinearMasterSlaveConstraint>(rConstraintId, master_dofs, slave_dofs, rRelationMatrix, rTranslationVector);
    rConstraintId++;
    return pconstraint;
}

/***********************************************************************************/
/***********************************************************************************/

void RVEPeriodicityUtility::Finalize(const Variable<array_1d<double, 3>>& rVariable)
{
    // First get the problem domain size
    const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Get the components
    std::vector<std::string> comp_list ({"_X","_Y"});
    if (domain_size == 3) {
        comp_list.push_back("_Z");
    }
    const std::string& r_base_variable_name = rVariable.Name();
    std::vector<const Variable<double>*> var_comp_vect(domain_size);
    for (std::size_t i = 0; i < domain_size; ++i) {
        var_comp_vect[i] = &KratosComponents<Variable<double>>::Get(r_base_variable_name + comp_list[i]);
    }

    for (auto& r_data : mAuxPairings) {
        auto& r_master_data = r_data.second;
        auto& r_master_ids = std::get<0>(r_master_data);
        auto& r_master_weights = std::get<1>(r_master_data);
        auto& r_T = std::get<2>(r_master_data);

        std::vector<IndexType> final_master_ids;
        std::vector<double> final_master_weights;
        Vector final_T = r_T;

        for (IndexType i = 0; i < r_master_ids.size(); ++i) {
            AppendIdsAndWeights(mAuxPairings, r_master_ids[i], r_master_weights[i], final_master_ids, final_master_weights, final_T);
        }

        // Assign back the finalized pairings and weights to the data structure
        r_master_ids = final_master_ids;
        r_master_weights = final_master_weights;
        r_T = final_T;
    }

    // First assign master and slave all to false
    auto& r_nodes_array = mrModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
        auto it_node = it_node_begin + i_node;
        it_node->Set(SLAVE, false);
        it_node->Set(MASTER, false);
    }
    // Compute the max id of the constraint
    IndexType constraint_id = 0;
    if (mrModelPart.NumberOfMasterSlaveConstraints() != 0) {
        constraint_id = (mrModelPart.MasterSlaveConstraints().end() - 1)->Id();
    }
    constraint_id++;

    // Define translation vector
    Vector aux_translation(1);

    ModelPart::MasterSlaveConstraintContainerType constraints;

    for (const auto& r_data : mAuxPairings) {
        const IndexType slave_id = r_data.first;
        const auto& r_master_data = r_data.second;
        auto& r_master_ids = std::get<0>(r_master_data);
        auto& r_master_weights = std::get<1>(r_master_data);
        auto& r_T = std::get<2>(r_master_data);

        // Very useful for debugging
        if (mEchoLevel > 0) {
            std::cout << "slave_id "  << slave_id << " - " << "master_ids ";
            for(auto& master_id : r_master_ids)
                std::cout << master_id << " ";
            std::cout << " - " << "master_weights ";
            for(auto& w : r_master_weights)
                std::cout << w << " " << "T ";
            std::cout << " - " << r_T << std::endl;
        }

        // Flag slave and master nodes
        mrModelPart.pGetNode(slave_id)->Set(SLAVE);
        for (auto id : r_master_ids) {
            mrModelPart.pGetNode(id)->Set(MASTER);
        }

        // Obtain the slave node
        auto pslave_node = mrModelPart.pGetNode(slave_id);

        // Define relation matrix (same for the different components)
        Matrix relation_matrix(1, r_master_weights.size());
        for (IndexType i = 0; i < relation_matrix.size2(); ++i) {
            relation_matrix(0, i) = r_master_weights[i];
        }

        // Generate the vector components constraints
        for (std::size_t i = 0; i < domain_size; ++i) {
            aux_translation[0] = r_T[i];
            constraints.push_back(GenerateConstraint(constraint_id, *var_comp_vect[i], pslave_node, r_master_ids, relation_matrix, aux_translation));
        }
    }

    mrModelPart.AddMasterSlaveConstraints(constraints.begin(), constraints.end());
}

}  // namespace Kratos.
