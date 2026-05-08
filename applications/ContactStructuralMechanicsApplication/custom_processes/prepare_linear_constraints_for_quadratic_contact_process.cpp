// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "custom_processes/prepare_linear_constraints_for_quadratic_contact_process.h"

namespace Kratos
{

void PrepareLinearConstraintsForQuadraticContactProcess::ExecuteInitialize()
{
    int max_condition_id = 0;
    auto& r_parent_mdpa = mrModelPart.GetParentModelPart();
    for (const auto& r_cond : r_parent_mdpa.MasterSlaveConstraints()) {
        if (r_cond.Id() > max_condition_id)
            max_condition_id = r_cond.Id();
    }
    max_condition_id++;

    for (auto& r_node : mrModelPart.Nodes()) {
        r_node.Set(VISITED, false);
    }

    std::vector<std::vector<std::size_t>> master_slave_nodes_ids = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
    MasterSlaveConstraint::DofPointerVectorType master_dofs(2);
    MasterSlaveConstraint::DofPointerVectorType slave_dofs(1);

    Matrix relation_matrix(1, 2);
    relation_matrix(0, 0) = 0.5;
    relation_matrix(0, 1) = 0.5;
    Vector constant_vector(2);
    constant_vector[0] = 0.0;
    constant_vector[1] = 0.0;

    for (auto& r_geom : mrModelPart.Geometries()) {

        IndexType count = 0;
        const IndexType number_of_points = r_geom.PointsNumber();
        for (IndexType it_node = 0; it_node < number_of_points; ++it_node) {
            if (r_geom[it_node].IsDefined(INTERFACE))
                if (r_geom[it_node].Is(INTERFACE))
                    ++count;
        }

        if (count == 6) { // quadratic triangles
            for (IndexType edge = 0; edge < 3; ++edge) {
                if (r_geom[master_slave_nodes_ids[edge][2]].IsNot(VISITED)) {

                    master_dofs[0] = r_geom[master_slave_nodes_ids[edge][0]].pGetDof(DISPLACEMENT_X);
                    master_dofs[1] = r_geom[master_slave_nodes_ids[edge][1]].pGetDof(DISPLACEMENT_X);
                    slave_dofs[0]  = r_geom[master_slave_nodes_ids[edge][2]].pGetDof(DISPLACEMENT_X);
                    mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",
                        max_condition_id, master_dofs, slave_dofs, relation_matrix, constant_vector);
                    ++max_condition_id;

                    master_dofs[0] = r_geom[master_slave_nodes_ids[edge][0]].pGetDof(DISPLACEMENT_Y);
                    master_dofs[1] = r_geom[master_slave_nodes_ids[edge][1]].pGetDof(DISPLACEMENT_Y);
                    slave_dofs[0]  = r_geom[master_slave_nodes_ids[edge][2]].pGetDof(DISPLACEMENT_Y);
                    mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",
                        max_condition_id, master_dofs, slave_dofs, relation_matrix, constant_vector);
                    ++max_condition_id;

                    master_dofs[0] = r_geom[master_slave_nodes_ids[edge][0]].pGetDof(DISPLACEMENT_Z);
                    master_dofs[1] = r_geom[master_slave_nodes_ids[edge][1]].pGetDof(DISPLACEMENT_Z);
                    slave_dofs[0] = r_geom[master_slave_nodes_ids[edge][2]].pGetDof(DISPLACEMENT_Z);
                    mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",
                        max_condition_id, master_dofs, slave_dofs, relation_matrix, constant_vector);
                    ++max_condition_id;

                    r_geom[master_slave_nodes_ids[edge][2]].Set(VISITED, true);
                }
            }
        }
    }
    // for (auto& r_node : mrModelPart.Nodes()) {
    //     r_node.Set(VISITED, false);
    // }
}

} // namespace Kratos
