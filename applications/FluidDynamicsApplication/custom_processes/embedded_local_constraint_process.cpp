//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes
#include <cstddef>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "embedded_local_constraint_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/

EmbeddedLocalConstraintProcess::EmbeddedLocalConstraintProcess(
    Model& rModel,
    Parameters& rParameters)
    : Process()
{
    // Validate input settings with defaults
    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Retrieve the required model parts
    const std::string model_part_name = rParameters["model_part_name"].GetString();
    mpModelPart = &rModel.GetModelPart(model_part_name).GetSubModelPart("fluid_computational_model_part");

    // Set whether constraints need to be adapted at each time step (FSI)
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mConstraintsAreCalculated = false;

    // Set to which elements constraints are applied
    mApplyToAllNegativeCutNodes = rParameters["apply_to_all_negative_cut_nodes"].GetBool();

    // Set whether intersection points should be added to the cloud points or only positive nodes
    mIncludeIntersectionPoints = rParameters["include_intersection_points"].GetBool();
    mSlipLength = rParameters["slip_length"].GetDouble();

    // Set whether nodal distances will be modified to avoid level set zeros
    mAvoidZeroDistances = rParameters["avoid_zero_distances"].GetBool();

    // Set which elements will not be active
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();

    // Define dofs (components) of a node and variables to constrain
    const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
    switch (n_dim) {
        case 2:
            mComponents = {"VELOCITY_X","VELOCITY_Y","PRESSURE"};
            mVariables = {&VELOCITY_X, &VELOCITY_Y};
            break;
        case 3:
            mComponents = {"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"};
            mVariables = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
            break;
        default:
            KRATOS_ERROR << "Wrong domain size.";
    }
}

void EmbeddedLocalConstraintProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void EmbeddedLocalConstraintProcess::ExecuteInitialize()
{
    KRATOS_TRY;

    // Continuous distance modification historical variables check
    KRATOS_ERROR_IF_NOT(mpModelPart->HasNodalSolutionStepVariable(DISTANCE)) << "DISTANCE variable not found in solution step variables list in " << mpModelPart->FullName() << ".\n";

    //TODO check for nodal variables (APPLY_EMBEDDED_CONSTRAINTS, EMBEDDED_CONSTRAINT_MASTERS, EMBEDDED_CONSTRAINT_MASTER_WEIGHTS)??

    KRATOS_CATCH("");
}

void EmbeddedLocalConstraintProcess::ExecuteBeforeSolutionLoop()
{
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void EmbeddedLocalConstraintProcess::ExecuteInitializeSolutionStep()
{
    if(!mConstraintsAreCalculated){
        // Avoid level set zero distances
        if (mAvoidZeroDistances) { ModifyDistances(); }

        // Deactivate full negative elements and fix their nodal values to zero
        if (mNegElemDeactivation) { this->DeactivateFullNegativeElements(); }

        // Declare support point weights and offsets for negative nodes of cut elements
        NodesCloudMapType clouds_map;
        NodesOffsetMapType offsets_map;

        // Calculate the negative nodes' clouds (support points)
        CalculateNodeClouds(clouds_map, offsets_map);

        // Hand over master node pointers and weights to slave nodes (for element)
        SetNodalValues(clouds_map);

        // Fix the slave dofs and calculate their values based on their master weights and offset (Dirichlet BC)
        FixAndCalculateConstrainedDofs(clouds_map, offsets_map);

        //OLD Apply a constraint to negative nodes of split elements on system level
        //ApplyConstraints(clouds_map, offsets_map);
        mConstraintsAreCalculated = true;
    }
}

void EmbeddedLocalConstraintProcess::ExecuteFinalizeSolutionStep()
{
    // Restore the initial state if the distance is checked at each time step
    if (mCheckAtEachStep){
        // Restore the flag to false
        mConstraintsAreCalculated = false;
        if (mNegElemDeactivation){
            // Recover the state previous to the element deactivation
            this->RecoverDeactivationPreviousState();
        }
        // Free constrained dofs
        for (const auto var : mVariables) {
            block_for_each(mSlaveNodes,[&](NodeType::Pointer p_slave_node){
                p_slave_node->Free(*var);
            });
        }
        mSlaveNodes.clear();
    }
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

void EmbeddedLocalConstraintProcess::CalculateNodeClouds(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    //Pointer for which function to call
    void(Kratos::EmbeddedLocalConstraintProcess::*fptrNodeClouds)(NodesCloudMapType&, NodesOffsetMapType&, std::vector<NodeType::Pointer>, std::vector<NodeType::Pointer>);
    if (mIncludeIntersectionPoints) {
        fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddAveragedNodeCloudsIncludingBC;
    } else {
        fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddAveragedNodeClouds;
    }

    // Containers for negative and positive side nodes of split element
    std::vector<NodePointerVectorType> neg_nodes;
    std::vector<NodePointerVectorType> pos_nodes;

    // Loop through all elements to get negative and positive nodes of split elements
    for (auto& rElement : mpModelPart->Elements()) {
    auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            if (mApplyToAllNegativeCutNodes || IsSmallCut(r_geom)) {
                // Containers for negative and positive side nodes of split element
                NodePointerVectorType neg_nodes_element;
                NodePointerVectorType pos_nodes_element;
                for (auto& rNode : r_geom) {
                    if (rNode.FastGetSolutionStepValue(DISTANCE) > 0.0) {
                        pos_nodes_element.push_back(&rNode);
                    } else {
                        neg_nodes_element.push_back(&rNode);
                    }
                }
                neg_nodes.push_back(neg_nodes_element);
                pos_nodes.push_back(pos_nodes_element);
            }
        }
    }

    mNSmallCutPos = 0;
    std::vector<NodePointerVectorType> neg_nodes_temp;
    std::vector<NodePointerVectorType> pos_nodes_temp;

    // Loop through cut elements to add single negative nodes first, then for two negative nodes of cut elements and so on
    std::size_t n_nodes = mpModelPart->Elements().begin()->GetGeometry().size();
    for (std::size_t n_neg_nodes = 1; n_neg_nodes < n_nodes; ++n_neg_nodes) {
        for(std::size_t i_elem = 0; i_elem < neg_nodes.size(); ++i_elem) {
            if (neg_nodes[i_elem].size() == n_neg_nodes) {
                (this->*fptrNodeClouds)(rCloudsMap, rOffsetsMap, neg_nodes[i_elem], pos_nodes[i_elem]);
            } else {
                neg_nodes_temp.push_back(neg_nodes[i_elem]);
                pos_nodes_temp.push_back(pos_nodes[i_elem]);
            }
        }
        neg_nodes = neg_nodes_temp;
        pos_nodes = pos_nodes_temp;
        neg_nodes_temp.clear();
        pos_nodes_temp.clear();
    }
    KRATOS_WATCH(mNSmallCutPos);
}

void EmbeddedLocalConstraintProcess::AddAveragedNodeClouds(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap,
    std::vector<NodeType::Pointer> neg_nodes_element,
    std::vector<NodeType::Pointer> pos_nodes_element)
{
    for (auto p_slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(p_slave_node)) {
            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();
            // Calculate weight as averaging of all master nodes
            const double cloud_node_weight = 1.0 / n_cloud_nodes;

            // Save positive element nodes (cloud nodes) and weights (averaged)
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], cloud_node_weight);
                cloud_data_vector(i_pos_node) = i_data;
            }

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(p_slave_node, 0.0));
            rCloudsMap.insert(std::make_pair(p_slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local averaged cloud for a negative node of a cut element");
        }
    }
}

void EmbeddedLocalConstraintProcess::AddAveragedNodeCloudsIncludingBC(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap,
    std::vector<NodeType::Pointer> neg_nodes_element,
    std::vector<NodeType::Pointer> pos_nodes_element)
{
    // Scale enforcement of boundary condition using the slip length
    const double value_gamma = 0.0;
    const double slip_length_mod = (mSlipLength < 1.0) ? std::abs(mSlipLength) : 1.0;

    for (auto p_slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(p_slave_node)) {
            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();

            // Get nodal distances
            double dist_slave = p_slave_node->FastGetSolutionStepValue(DISTANCE);
            double sum_dist_cloud_nodes = 0.0;
            for (auto pos_node : pos_nodes_element) {
                sum_dist_cloud_nodes += pos_node->FastGetSolutionStepValue(DISTANCE);
            }
            //TODO: case: small cut for positive side - change!!!  // TODO: scale slave_distance using edge length/elements size?
            // --> apply boundary condition by making negative distance very small
            // TODO: tolerance ??? -- 1e-6 IsSmallCut, 1e-9 sufficient (?), 1e-10 too small!
            if (sum_dist_cloud_nodes < 1e-8) {
                dist_slave = sum_dist_cloud_nodes;
                sum_dist_cloud_nodes = n_cloud_nodes;
                mNSmallCutPos++;
            }

            // Calculate weight
            //without slip length:
            //const double cloud_node_weight = dist_slave / sum_dist_cloud_nodes;
            const double cloud_node_weight = slip_length_mod / n_cloud_nodes + dist_slave / sum_dist_cloud_nodes * (1.0 - slip_length_mod);
            //KRATOS_WATCH(cloud_node_weight);

            // Save positive element nodes (cloud nodes) and weights (linear interpolation of averaged positive nodes and averaged intersection points)
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], cloud_node_weight);
                cloud_data_vector(i_pos_node) = i_data;
            }

            //Calculate offset for master-slave constraint
            //without slip length:
            //const double offset = value_gamma * (1 - dist_slave * n_cloud_nodes / sum_dist_cloud_nodes);
            const double offset = value_gamma * (1.0 - slip_length_mod) * (1.0 - dist_slave * n_cloud_nodes / sum_dist_cloud_nodes);

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(p_slave_node, offset));
            rCloudsMap.insert(std::make_pair(p_slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local averaged cloud including intersection points for a negative node of a cut element");
        }
    }
}

void EmbeddedLocalConstraintProcess::SetNodalValues(NodesCloudMapType& rCloudsMap)
{
    // Set APPLY_EMBEDDED_CONSTRAINTS to false for all nodes
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.SetValue(APPLY_EMBEDDED_CONSTRAINTS, false);
    });

    const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];

    // Loop through all slave nodes
    // TODO run parallel?
    for (const auto& r_slave : rCloudsMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;

        // Create master vector and relation matrix for the negative node
        NodePointerVectorType master_nodes;
        MatrixType dof_relation_matrix;
        master_nodes.reserve(r_ext_op_data.size());
        dof_relation_matrix.resize(n_dim, r_ext_op_data.size()*n_dim, false);

        // Fill vector of pointers to master nodes and Matrix of master dof weights
        std::size_t i_master = 0;
        for (auto& r_cloud_data : r_ext_op_data) {
            const auto p_cloud_node = r_cloud_data.first;
            const double cloud_node_weight = r_cloud_data.second;

            // Fill vector of pointers to master nodes
            master_nodes.push_back(p_cloud_node);

            // Get weights for all dofs of the master node to build relation matrix of slave node
            for (std::size_t d1 = 0; d1 < n_dim; ++d1) {
                for (std::size_t d2 = 0; d2 < n_dim; ++d2) {
                    dof_relation_matrix(d1, i_master*n_dim+d2) = (d1 == d2) ? cloud_node_weight : 0.0;
                }
            }
            i_master++;
        }

        // Set nodal variables
        //TODO Synchronize them?? //mrModelPart.GetCommunicator().AssembleNonHistoricalData(EMBEDDED_IS_ACTIVE);
        p_slave_node->SetValue(APPLY_EMBEDDED_CONSTRAINTS, true);
        p_slave_node->SetValue(EMBEDDED_CONSTRAINT_MASTERS, master_nodes);
        p_slave_node->SetValue(EMBEDDED_CONSTRAINT_MASTER_WEIGHTS, dof_relation_matrix);
    }
}

void EmbeddedLocalConstraintProcess::FixAndCalculateConstrainedDofs(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    // Loop through all slave nodes
    // TODO run parallel?
    for (const auto& r_slave : rCloudsMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;

        // Fix velocity dofs of constrained node and calculate their values
        for (const auto var : mVariables) {
            p_slave_node->Fix(*var);
            mSlaveNodes.push_back(p_slave_node);

            double dof_value = 0.0;
            for (auto& r_cloud_data : r_ext_op_data) {
                const auto p_cloud_node = r_cloud_data.first;
                const double cloud_node_weight = r_cloud_data.second;

                dof_value += cloud_node_weight * p_cloud_node->FastGetSolutionStepValue(*var);
            }
            p_slave_node->FastGetSolutionStepValue(*var) = dof_value;
        }
    }
}

void EmbeddedLocalConstraintProcess::ApplyConstraints(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    // Initialize counter of master slave constraints
    ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;

    // Loop through all negative nodes of split elements (slave nodes)
    for (const auto& r_slave : rCloudsMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;
        double offset = rOffsetsMap[p_slave_node];

        // Create vectors and relation matrix for the negative node
        DofPointerVectorType slave_dofs;
        DofPointerVectorType master_dofs;
        MatrixType dof_relation_matrix;
        VectorType dof_offset;

        const std::size_t n_var = mVariables.size();
        dof_offset.resize(n_var, false);
        dof_relation_matrix.resize(n_var, r_ext_op_data.size()*n_var, false);
        slave_dofs.reserve(n_var);
        master_dofs.reserve(r_ext_op_data.size()*n_var);

        // Get all dofs of the negative node and their constraint constants
        std::size_t i_var = 0;
        for (const auto var : mVariables) {
            slave_dofs.push_back(p_slave_node->pGetDof(*var));
            dof_offset(i_var++) = offset;
        }

        // Get all dofs of all cloud nodes of the negative node and their weights
        std::size_t i_cloud = 0;
        for (auto& r_cloud_data : r_ext_op_data) {
            const auto p_cloud_node = r_cloud_data.first;
            const double cloud_node_N = r_cloud_data.second;

            std::size_t i_cloud_var = 0;
            for (const auto var : mVariables) {
                master_dofs.push_back(p_cloud_node->pGetDof(*var));

                for (std::size_t i_var = 0; i_var < n_var; i_var++) {
                    dof_relation_matrix(i_var, i_cloud*n_var+i_cloud_var) = (i_var == i_cloud_var) ? cloud_node_N : 0.0;
                }
                i_cloud_var++;
            }
            i_cloud++;
        }

        // Create new linear master-slave constraint for the negative node
        // (It is faster to create constraints for all dofs of one node as one new constraint, instead of creating them all separately.)
        mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
        master_dofs, slave_dofs, dof_relation_matrix, dof_offset);
    }
}

void EmbeddedLocalConstraintProcess::DeactivateFullNegativeElements()
{
    ModelPart::NodesContainerType& rNodes = mpModelPart->Nodes();
    ModelPart::ElementsContainerType& rElements = mpModelPart->Elements();

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        ModelPart::NodesContainerType::iterator it_node = rNodes.begin() + i_node;
        it_node->SetValue(EMBEDDED_IS_ACTIVE, 0);
    }

    // Deactivate those elements whose negative distance nodes summation is equal to their number of nodes
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(rElements.size()); ++i_elem){
        unsigned int n_neg = 0;
        ModelPart::ElementsContainerType::iterator it_elem = rElements.begin() + i_elem;
        auto& rGeometry = it_elem->GetGeometry();

        // Check the distance function sign at the element nodes
        for (unsigned int i_node = 0; i_node<rGeometry.size(); i_node++){
            if (rGeometry[i_node].FastGetSolutionStepValue(DISTANCE) < 0.0){
                n_neg++;
            }
        }

        (n_neg == rGeometry.size()) ? it_elem->Set(ACTIVE, false) : it_elem->Set(ACTIVE, true);

        // If the element is ACTIVE, all its nodes are active as well
        if (it_elem->Is(ACTIVE)){
            for (unsigned int i_node = 0; i_node < rGeometry.size(); ++i_node){
                int& activation_index = rGeometry[i_node].GetValue(EMBEDDED_IS_ACTIVE);
                #pragma omp atomic
                activation_index += 1;
            }
        }
    }

    // Synchronize the EMBEDDED_IS_ACTIVE variable flag
    mpModelPart->GetCommunicator().AssembleNonHistoricalData(EMBEDDED_IS_ACTIVE);

    // Set to zero and fix the DOFs in the remaining inactive nodes
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        auto it_node = rNodes.begin() + i_node;
        if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
            for (std::size_t i_var = 0; i_var < mComponents.size(); i_var++){
                const auto& r_double_var = KratosComponents<Variable<double>>::Get(mComponents[i_var]);
                // Fix the nodal DOFs
                it_node->Fix(r_double_var);
                // Set to zero the nodal DOFs
                it_node->FastGetSolutionStepValue(r_double_var) = 0.0;
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::RecoverDeactivationPreviousState()
{
    // Activate again all the elements
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mpModelPart->NumberOfElements()); ++i_elem){
        auto it_elem = mpModelPart->ElementsBegin() + i_elem;
        it_elem->Set(ACTIVE, true);
    }
    // Free the negative DOFs that were fixed
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(mpModelPart->NumberOfNodes()); ++i_node){
        auto it_node = mpModelPart->NodesBegin() + i_node;
        if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
            for (std::size_t i_var = 0; i_var < mComponents.size(); i_var++){
                const auto& r_double_var = KratosComponents<Variable<double>>::Get(mComponents[i_var]);
                // Free the nodal DOFs  that were fixed
                it_node->Free(r_double_var);
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::ModifyDistances()
{
    auto& r_nodes = mpModelPart->Nodes();
    const double tol_d = 1.0e-11;

    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(r_nodes.size()); ++i_node) {
        auto it_node = r_nodes.begin() + i_node;
        double& d = it_node->FastGetSolutionStepValue(DISTANCE);

        // Check if the distance values are close to zero, if so set the tolerance as distance value
        if (std::abs(d) < tol_d) {
            d = (d > 0.0) ? tol_d : -tol_d;
        }
    }
}

bool EmbeddedLocalConstraintProcess::IsSplit(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    std::size_t n_pos = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        } else {
            n_pos++;
        }
    }
    return (n_pos != 0 && n_neg != 0);
}

bool EmbeddedLocalConstraintProcess::IsSmallCut(const GeometryType& rGeometry)
{
    const double tol_d = 0.000001;
    for (const auto& r_node : rGeometry) {
        if (std::abs(r_node.FastGetSolutionStepValue(DISTANCE)) < tol_d) {
            return true;
        }
    }
    return false;
}

bool EmbeddedLocalConstraintProcess::IsNegative(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        }
    }
    return (n_neg == rGeometry.PointsNumber());
}

}; // namespace Kratos.
