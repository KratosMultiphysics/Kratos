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
#include "rve_periodicity_utility.h"
#include "utilities/binbased_fast_point_locator_conditions.h"
#include "constraints/linear_master_slave_constraint.h"

namespace Kratos
{
void RVEPeriodicityUtility::AssignPeriodicity(
    ModelPart& rMasterModelPart,
    ModelPart& rSlaveModelPart,
    const Matrix& rStrainTensor,
    const Vector& rDirection
    )
{
    KRATOS_ERROR_IF(rMasterModelPart.NumberOfConditions() == 0) << "the master is expected to have conditions and it is empty" << std::endl;

    const Vector translation = prod(rStrainTensor, rDirection);

    BinBasedFastPointLocatorConditions<3> bin_based_point_locator(rMasterModelPart);
    bin_based_point_locator.UpdateSearchDatabase();

    int max_search_results = 100;
    double search_tolerance = 1e-6;

    // Construct auxiliary data structure to contain the master slave relation.
    // Slave nodes must appear once, however a non-circular dependency is allowed between the masters
    for (IndexType i = 0; i < rSlaveModelPart.Nodes().size(); ++i) {
        // Search in which condition it falls
        auto it_node = rSlaveModelPart.NodesBegin() + i;

        Condition::Pointer p_host_cond;
        Vector N;
        array_1d<double, 3> transformed_slave_coordinates = it_node->Coordinates() - rDirection;

        // Finding the host element for this node
        const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(transformed_slave_coordinates, N, p_host_cond, max_search_results, search_tolerance);
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
            KRATOS_ERROR << "Counterpart not found for slave node " << it_node->Id() << std::endl;
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
    // Get the components
    const std::string& r_base_variable_name = rVariable.Name();
    auto& r_var_x = KratosComponents<Variable<double>>::Get(r_base_variable_name + "_X");
    auto& r_var_y = KratosComponents<Variable<double>>::Get(r_base_variable_name + "_Y");
    auto& r_var_z = KratosComponents<Variable<double>>::Get(r_base_variable_name + "_Z");

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
    Vector xtranslation(1);
    Vector ytranslation(1);
    Vector ztranslation(1);

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

        xtranslation[0] = r_T[0];
        constraints.push_back(GenerateConstraint(constraint_id, r_var_x, pslave_node, r_master_ids, relation_matrix, xtranslation));

        ytranslation[0] = r_T[1];
        constraints.push_back(GenerateConstraint(constraint_id, r_var_y, pslave_node, r_master_ids, relation_matrix, ytranslation));

        ztranslation[0] = r_T[2];
        constraints.push_back(GenerateConstraint(constraint_id, r_var_z, pslave_node, r_master_ids, relation_matrix, ztranslation));
    }

    mrModelPart.AddMasterSlaveConstraints(constraints.begin(), constraints.end());
}

}  // namespace Kratos.


