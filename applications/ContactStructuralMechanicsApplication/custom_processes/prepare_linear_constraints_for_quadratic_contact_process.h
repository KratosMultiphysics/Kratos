// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Your Name
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PrepareLinearConstraintsForQuadraticContactProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Description of the process for preparing linear constraints for quadratic contact
 * @author Your Name
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PrepareLinearConstraintsForQuadraticContactProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// General type definitions
    using NodesArrayType = ModelPart::NodesContainerType;

    /// Pointer definition of PrepareLinearConstraintsForQuadraticContactProcess
    KRATOS_CLASS_POINTER_DEFINITION( PrepareLinearConstraintsForQuadraticContactProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the process
     * @param rModelPart The model part to be considered
     * @param rProcessInfo The process info to be considered
     */
    PrepareLinearConstraintsForQuadraticContactProcess(
        ModelPart& rModelPart,
        ProcessInfo& rProcessInfo
        ) : mrModelPart(rModelPart),
            mrProcessInfo(rProcessInfo)
    {
    }

    virtual ~PrepareLinearConstraintsForQuadraticContactProcess() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms
     */
    void Execute() override
    {
    }

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override
    {
        int max_condition_id = 1000; // TODO take the actual max id from the model part

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
        constant_vector.clear();

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
                        mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", max_condition_id, master_dofs,
                            slave_dofs, relation_matrix, constant_vector);
                        ++max_condition_id;
                        master_dofs[0] = r_geom[master_slave_nodes_ids[edge][0]].pGetDof(DISPLACEMENT_Y);
                        master_dofs[1] = r_geom[master_slave_nodes_ids[edge][1]].pGetDof(DISPLACEMENT_Y);
                        slave_dofs[0]  = r_geom[master_slave_nodes_ids[edge][2]].pGetDof(DISPLACEMENT_Y);
                        mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", max_condition_id, master_dofs,
                            slave_dofs, relation_matrix, constant_vector);
                        ++max_condition_id;
    
                        master_dofs[0] = r_geom[master_slave_nodes_ids[edge][0]].pGetDof(DISPLACEMENT_Z);
                        master_dofs[1] = r_geom[master_slave_nodes_ids[edge][1]].pGetDof(DISPLACEMENT_Z);
                        slave_dofs[0] = r_geom[master_slave_nodes_ids[edge][2]].pGetDof(DISPLACEMENT_Z);
                        mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", max_condition_id, master_dofs,
                            slave_dofs, relation_matrix, constant_vector);
                        ++max_condition_id;

                        r_geom[master_slave_nodes_ids[edge][2]].Set(VISITED, true);
                    }
                }
            }
        }
    }

    /**
     * @brief This function is designed for being execute once before the solution loop but after
     * all of the solvers where built
     */
    void ExecuteBeforeSolutionLoop() override
    {
    }

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
    }

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
    }

    /**
     * @brief This function will be executed at every time step BEFORE  writing the output
     */
    void ExecuteBeforeOutputStep() override
    {
    }

    /**
     * @brief This function will be executed at every time step AFTER writing the output
     */
    void ExecuteAfterOutputStep() override
    {
    }

    /**
     * @brief This function is designed for being called at the end of the computations
     */
    void ExecuteFinalize() override
    {
    }

    /**
     * @brief This function is designed for being called after ExecuteInitialize ONCE
     * to verify that the input is correct.
     */
    int Check() override
    {
        return 0;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        return Parameters(R"({})");
    }

    ///@}

private:
    ///@name Name and ID
    ///@{

    ///@}
    ///@name Data attributes
    ///@{

    /// The model part reference
    ModelPart& mrModelPart;

    /// The process info reference
    ProcessInfo& mrProcessInfo;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Private access
    ///@{

    ///@}
    ///@name Unallowed copy and assignment
    ///@{

    /// Copy constructor
    PrepareLinearConstraintsForQuadraticContactProcess(const PrepareLinearConstraintsForQuadraticContactProcess&) = default;

    /// Assignment operator
    PrepareLinearConstraintsForQuadraticContactProcess& operator=(const PrepareLinearConstraintsForQuadraticContactProcess&) = delete;

    ///@}

}; /// Class PrepareLinearConstraintsForQuadraticContactProcess.

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
///@name Search methods
///@{

///@}
///@name Friends
///@{

///@}

///@}
///@name Friends
///@{

///@}

} // namespace Kratos
