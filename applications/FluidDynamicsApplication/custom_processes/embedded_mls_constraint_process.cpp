//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Franziska Wahl
//

// System includes
#include <cmath>
#include <unordered_set>

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "embedded_mls_constraint_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

/* Public functions *******************************************************/

EmbeddedMLSConstraintProcess::EmbeddedMLSConstraintProcess(
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
    mCheckAtEachStep = rParameters["check_at_each_step"].GetBool();
    mConstraintsAreCalculated = false;

    // Set to which elements constraints are applied
    mApplyToAllNegativeCutNodes = rParameters["apply_to_all_negative_cut_nodes"].GetBool();

    // Set the order of the MLS extension operator used in the MLS shape functions utility
    mMLSExtensionOperatorOrder = rParameters["mls_extension_operator_order"].GetInt();

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

void EmbeddedMLSConstraintProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void EmbeddedMLSConstraintProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    // Continuous distance modification historical variables check
    KRATOS_ERROR_IF_NOT(mpModelPart->HasNodalSolutionStepVariable(DISTANCE)) << "DISTANCE variable not found in solution step variables list in " << mpModelPart->FullName() << ".\n";
    KRATOS_CATCH("");
}

void EmbeddedMLSConstraintProcess::ExecuteBeforeSolutionLoop()
{
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void EmbeddedMLSConstraintProcess::ExecuteInitializeSolutionStep()
{
    if(!mConstraintsAreCalculated){
        // Avoid level set zero distances
        if (mAvoidZeroDistances) { ModifyDistances(); }

        // Set the required interface flags
        // NOTE: this will deactivate all negative as well as intersected elements
        SetInterfaceFlags();

        // Declare support point weights and offsets for negative nodes of cut elements
        NodesCloudMapType clouds_map;
        NodesOffsetMapType offsets_map;

        // Calculate the negative nodes' clouds (support points) as extension basis
        if (mIncludeIntersectionPoints) {
            CalculateNodeCloudsIncludingBC(clouds_map, offsets_map);
        } else {
            CalculateNodeClouds(clouds_map, offsets_map);
        }

        // Reactivate intersected and negative elements depending on process settings
        // Fix nodal values of deactivated nodes to zero
        ReactivateElementsAndFixNodes();

        // Hand over master node pointers and weights to slave nodes (for element)
        SetNodalValues(clouds_map);

        // Fix the slave dofs and calculate their values based on their master weights and offset (Dirichlet BC)
        FixAndCalculateConstrainedDofs(clouds_map, offsets_map);

        // Store master nodes, weights and offsets for subsequent iterations
        mCloudsMap = clouds_map;
        mOffsetsMap = offsets_map;
        mConstraintsAreCalculated = true;
    } else {
        UpdateConstrainedDofs();
    }
}

void EmbeddedMLSConstraintProcess::ExecuteFinalizeSolutionStep()
{
    // Restore the initial state if the distance is checked at each time step
    if (mCheckAtEachStep){
        // Restore the flag to false
        mConstraintsAreCalculated = false;
        if (mNegElemDeactivation){
            // Recover the state previous to the element deactivation
            this->RecoverDeactivationPreviousState();
        }

        FreeConstrainedDofs();

        mCloudsMap.clear();
        mOffsetsMap.clear();
    }
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

void EmbeddedMLSConstraintProcess::CalculateNodeClouds(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    // Get the MLS shape functions function
    auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

    // Loop the elements to create the negative nodes MLS basis
    for (auto& rElement : mpModelPart->Elements()) {
        // Check if the element is split
        const auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            if (mApplyToAllNegativeCutNodes || IsSmallCut(r_geom)) {
                // Find the intersected element negative nodes
                for (auto& r_node : r_geom) {
                    // Check whether node is negative (not ACTIVE) or already was added because of a previous element
                    const std::size_t found = rCloudsMap.count(&r_node);
                    if (r_node.IsNot(ACTIVE) && !found) {

                        // Get the current negative node neighbours cloud
                        Matrix cloud_nodes_coordinates;
                        PointerVector<NodeType> cloud_nodes;
                        PointerVector<NodeType> positive_neighbor_nodes;
                        SetNegativeNodeSupportCloud(r_node, cloud_nodes, positive_neighbor_nodes, cloud_nodes_coordinates);

                        // Calculate the MLS basis in the current negative node
                        Vector N_container;
                        const array_1d<double,3> r_coords = r_node.Coordinates();
                        const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                        p_mls_sh_func(cloud_nodes_coordinates, r_coords, 1.01 * mls_kernel_rad, N_container);

                        // Save the extension operator nodal data
                        std::size_t n_cl_nod = cloud_nodes.size();
                        CloudDataVectorType cloud_data_vector(n_cl_nod);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cl_nod; ++i_cl_nod) {
                            auto p_cl_node = cloud_nodes(i_cl_nod);
                            auto i_data = std::make_pair(p_cl_node, N_container[i_cl_nod]);
                            cloud_data_vector(i_cl_nod) = i_data;
                        }

                        // Pair pointer to negative node with vector of pairs of a pointer to a positive node with its MLS shape function value and add it to the extension operator map
                        auto ext_op_key_data = std::make_pair(&r_node, cloud_data_vector);
                        rCloudsMap.insert(ext_op_key_data);

                        // Pair pointer to negative node with offset for the master-slave constraint of that node and add it to the constraint offset map
                        auto offset_key_data = std::make_pair(&r_node, 0.0);
                        rOffsetsMap.insert(offset_key_data);

                        KRATOS_WATCH("Added mls cloud for a negative node of a cut element");
                    }
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::CalculateNodeCloudsIncludingBC(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    // Get the MLS shape functions function
    auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

    // Scale enforcement of boundary condition using the slip length
    const double value_gamma = 0.0;
    const double slip_length_mod = mSlipLength; //(mSlipLength < 1.0) ? std::abs(mSlipLength) : 1.0;

    // Loop the elements to create the negative nodes MLS basis
    for (auto& rElement : mpModelPart->Elements()) {
        // Check if the element is split
        const auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            if (mApplyToAllNegativeCutNodes || IsSmallCut(r_geom)) {
                // Find the intersected element negative nodes
                for (auto& r_node : r_geom) {
                    // Check whether node is negative (not ACTIVE) or already was added because of a previous element
                    const std::size_t found = rCloudsMap.count(&r_node);
                    if (r_node.IsNot(ACTIVE) && !found) {

                        // Get the current negative node neighbors cloud and directly neighboring positive nodes
                        Matrix cloud_nodes_coordinates;
                        PointerVector<NodeType> cloud_nodes;
                        PointerVector<NodeType> positive_neighbor_nodes;
                        SetNegativeNodeSupportCloud(r_node, cloud_nodes, positive_neighbor_nodes, cloud_nodes_coordinates);

                        // Get the current negative node intersection points
                        Matrix intersection_points_coordinates;
                        SetNegativeNodeIntersectionPoints(r_node, positive_neighbor_nodes, intersection_points_coordinates);

                        // Calculate Kernal Radius using the cloud nodes only (max. node distance to the negative node)
                        const array_1d<double,3> r_coords = r_node.Coordinates();
                        const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);

                        // Scale weighting of intersection points by modifying their distance to the negative node based on the slip length
                        const std::size_t n_int_pt = positive_neighbor_nodes.size();
                        for (std::size_t i_int_pt = 0; i_int_pt < n_int_pt; ++i_int_pt) {
                            // Calculate the vector pointing to the negative node from the intersection point and its length
                            array_1d<double,3> slip_vector;
                            slip_vector(0) = r_coords(0) - intersection_points_coordinates(i_int_pt, 0);
                            slip_vector(1) = r_coords(1) - intersection_points_coordinates(i_int_pt, 1);
                            slip_vector(2) = r_coords(2) - intersection_points_coordinates(i_int_pt, 2);
                            const double slip_vector_length = std::sqrt( std::pow(slip_vector(0),2) + std::pow(slip_vector(1),2) + std::pow(slip_vector(2),2) );
                            // Scale the vector using a distance to the negative node of twice the kernel radius as minimal full slip condition (slip_length=1.0)
                            slip_vector(0) = 2.0 * slip_vector(0) + slip_length_mod * (2.0 * mls_kernel_rad - slip_vector_length) / slip_vector_length;
                            slip_vector(1) = 2.0 * slip_vector(0) + slip_length_mod * (2.0 * mls_kernel_rad - slip_vector_length) / slip_vector_length;
                            slip_vector(2) = 2.0 * slip_vector(0) + slip_length_mod * (2.0 * mls_kernel_rad - slip_vector_length) / slip_vector_length;
                            // Modify the distance of the intersection point and thereby its weighting
                            intersection_points_coordinates(i_int_pt, 0) += slip_vector(0);
                            intersection_points_coordinates(i_int_pt, 1) += slip_vector(1);
                            intersection_points_coordinates(i_int_pt, 2) += slip_vector(2);
                        }

                        // Add matrix of cloud nodes to matrix of intersection points to get matrix of all cloud points
                        Matrix cloud_points_coordinates;
                        const std::size_t n_cl_nod = cloud_nodes.size();
                        cloud_points_coordinates.resize(n_cl_nod+n_int_pt, 3);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cl_nod; ++i_cl_nod) {
                            cloud_points_coordinates(i_cl_nod, 0) = cloud_nodes_coordinates(i_cl_nod, 0);
                            cloud_points_coordinates(i_cl_nod, 1) = cloud_nodes_coordinates(i_cl_nod, 1);
                            cloud_points_coordinates(i_cl_nod, 2) = cloud_nodes_coordinates(i_cl_nod, 2);
                        }
                        for (std::size_t i_int_pt = 0; i_int_pt < n_int_pt; ++i_int_pt) {
                            cloud_points_coordinates(n_cl_nod+i_int_pt, 0) = intersection_points_coordinates(i_int_pt, 0);
                            cloud_points_coordinates(n_cl_nod+i_int_pt, 1) = intersection_points_coordinates(i_int_pt, 1);
                            cloud_points_coordinates(n_cl_nod+i_int_pt, 2) = intersection_points_coordinates(i_int_pt, 2);
                        }

                        // Calculate the MLS basis in the current negative node
                        Vector N_container;
                        p_mls_sh_func(cloud_points_coordinates, r_coords, 1.01 * mls_kernel_rad, N_container);

                        // Save the extension operator nodal data
                        CloudDataVectorType cloud_data_vector(n_cl_nod);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cl_nod; ++i_cl_nod) {
                            auto p_cl_node = cloud_nodes(i_cl_nod);
                            auto i_data = std::make_pair(p_cl_node, N_container[i_cl_nod]);
                            cloud_data_vector(i_cl_nod) = i_data;
                        }

                        // Calculate Master-slave constraint offset based on BC values at the intersection points and their MLS shape function values
                        double constraint_offset = 0.0;
                        for (std::size_t i_int_pt = 0; i_int_pt < n_int_pt; ++i_int_pt) {
                            // Get boundary condition value for the intersection point  // TODO  different boundary condition?!?
                            //const double temp_bc = rElement.GetValue(EMBEDDED_SCALAR);  // TODO =0.0 if value_bc = 0.0 ???
                            //const double value_bc = std::pow(intersection_points_coordinates(i_int_pt, 0), 2) + std::pow(intersection_points_coordinates(i_int_pt, 1), 2);
                            // Add influence of value at intersection point to the constraint for the negative node
                            constraint_offset += value_gamma * N_container[n_cl_nod+i_int_pt];
                        }

                        // Pair pointer to negative node with vector of pairs of a pointer to a positive node with its MLS shape function value and add it to the extension operator map
                        auto ext_op_key_data = std::make_pair(&r_node, cloud_data_vector);
                        rCloudsMap.insert(ext_op_key_data);

                        // Pair pointer to negative node with offset for the master-slave constraint of that node and add it to the constraint offset map
                        auto offset_key_data = std::make_pair(&r_node, constraint_offset);
                        rOffsetsMap.insert(offset_key_data);

                        KRATOS_WATCH("Added mls cloud including intersection points for a negative node of a cut element");
                    }
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::SetNodalValues(NodesCloudMapType& rCloudsMap)
{
    // Set flag APPLY_EMBEDDED_CONSTRAINTS to false for all nodes
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.SetValue(APPLY_EMBEDDED_CONSTRAINTS, 0);
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
        p_slave_node->SetValue(APPLY_EMBEDDED_CONSTRAINTS, 1);
        p_slave_node->SetValue(EMBEDDED_CONSTRAINT_MASTERS, master_nodes);
        p_slave_node->SetValue(EMBEDDED_CONSTRAINT_MASTER_WEIGHTS, dof_relation_matrix);
    }
}

void EmbeddedMLSConstraintProcess::FixAndCalculateConstrainedDofs(
    NodesCloudMapType& rCloudsMap,
    NodesOffsetMapType& rOffsetsMap)
{
    // Loop through all slave nodes
    // TODO run parallel?
    for (const auto& r_slave : rCloudsMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;
        const double offset = rOffsetsMap[p_slave_node];

        // Fix velocity dofs of constrained node and calculate their values
        for (const auto var : mVariables) {
            p_slave_node->Fix(*var);

            double dof_value = offset;
            for (auto& r_cloud_data : r_ext_op_data) {
                const auto p_cloud_node = r_cloud_data.first;
                const double cloud_node_weight = r_cloud_data.second;

                dof_value += cloud_node_weight * p_cloud_node->FastGetSolutionStepValue(*var);
            }
            p_slave_node->FastGetSolutionStepValue(*var) = dof_value;
        }
    }
}

void EmbeddedMLSConstraintProcess::UpdateConstrainedDofs()
{
    // Loop through all slave nodes | TODO: run parallel?
    for (const auto& r_slave : mCloudsMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;
        const double offset = mOffsetsMap[p_slave_node];

        // Update the velocity values of the constrained node
        for (const auto var : mVariables) {
            double dof_value = offset;
            for (auto& r_cloud_data : r_ext_op_data) {
                const auto p_cloud_node = r_cloud_data.first;
                const double cloud_node_weight = r_cloud_data.second;
                dof_value += cloud_node_weight * p_cloud_node->FastGetSolutionStepValue(*var);
            }
            p_slave_node->FastGetSolutionStepValue(*var) = dof_value;
        }
    }
}

void EmbeddedMLSConstraintProcess::FreeConstrainedDofs()
{
    // Get all pointers to slave nodes as vector
    NodePointerVectorType slave_nodes(mCloudsMap.size());
    for (const auto& r_slave : mCloudsMap) {
        slave_nodes.push_back(r_slave.first);
    }

    // Free constrained dofs of slave nodes
    for (const auto var : mVariables) {
        block_for_each(slave_nodes,[&](NodeType::Pointer p_slave_node){
            p_slave_node->Free(*var);
        });
    }
}

void EmbeddedMLSConstraintProcess::SetInterfaceFlags()
{
    //TODO only set ACTIVE or not, no need for BOUNDARY and INTERFACE ?
    // Initialize flags to false
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.Set(ACTIVE, false); // Nodes that belong to the elements to be assembled
        rNode.Set(BOUNDARY, false); // Nodes that belong to the surrogate boundary
        rNode.Set(INTERFACE, false); // Nodes that belong to the cloud of points
    });
    block_for_each(mpModelPart->Elements(), [](Element& rElement){
        rElement.Set(ACTIVE, false); // Elements in the positive distance region (the ones to be assembled)
        rElement.Set(BOUNDARY, false); // Intersected elements
        rElement.Set(INTERFACE, false); // Positive distance elements owning the surrogate boundary nodes
    });

    // Find active regions
    for (auto& rElement : mpModelPart->Elements()) {
        auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            // Mark the intersected elements as BOUNDARY
            // The intersected elements are also flagged as non ACTIVE to avoid assembling them
            // Note that the split elements BC is applied by means of the extension operators
            rElement.Set(ACTIVE, false);
            rElement.Set(BOUNDARY, true);
            // Mark the positive intersected nodes as BOUNDARY
            for (auto& rNode : r_geom) {
                if (!(rNode.FastGetSolutionStepValue(DISTANCE) < 0.0)) {
                    rNode.Set(BOUNDARY, true);
                }
            }
        } else if (!IsNegative(r_geom)) {
            // Mark the positive element as ACTIVE to assemble it
            rElement.Set(ACTIVE, true);
            // Mark the active element nodes as ACTIVE
            for (auto& rNode : r_geom) {
                rNode.Set(ACTIVE, true);
            }
        }
    }

    // Find the surrogate boundary elements
    // Note that we rely on the fact that the neighbours are sorted according to the faces
    for (auto& rElement : mpModelPart->Elements()) {
        if (rElement.Is(BOUNDARY)) {
            const auto& r_geom = rElement.GetGeometry();
            const std::size_t n_faces = r_geom.FacesNumber();
            auto& r_neigh_elems = rElement.GetValue(NEIGHBOUR_ELEMENTS);
            for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
                // The neighbour corresponding to the current face is ACTIVE means that the current face is surrogate boundary
                // Flag the current neighbour owning the surrogate face as INTERFACE
                auto& r_neigh_elem = r_neigh_elems[i_face];
                if (r_neigh_elem.Is(ACTIVE)) {
                    r_neigh_elem.Set(INTERFACE, true);
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::ReactivateElementsAndFixNodes()
{
    // Make intersected elements and their nodes ACTIVE
    for (auto& rElement : mpModelPart->Elements()) {
        // Check if the element is split
        const auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            rElement.Set(ACTIVE, true);
            for (auto& rNode : r_geom) {
                rNode.Set(ACTIVE, true);
            }
        }
    }

    // Make negative elements and their nodes ACTIVE
    if ( !mNegElemDeactivation ) {
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is negative
            const auto& r_geom = rElement.GetGeometry();
            if (IsNegative(r_geom)) {
                rElement.Set(ACTIVE, true);
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, true);
                }
            }
        }
    } else {
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
        for (int k = 0; k < static_cast<int>(rElements.size()); ++k){
            ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
            auto& rGeometry = itElement->GetGeometry();

            // If the element is ACTIVE, all its nodes are active as well
            if (itElement->Is(ACTIVE)){
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
}

void EmbeddedMLSConstraintProcess::RecoverDeactivationPreviousState()
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

void EmbeddedMLSConstraintProcess::ModifyDistances()
{
    auto& r_nodes = mpModelPart->Nodes();
    double tol_d = 1.0e-11;

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

bool EmbeddedMLSConstraintProcess::IsSplit(const GeometryType& rGeometry)
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

bool EmbeddedMLSConstraintProcess::IsSmallCut(const GeometryType& rGeometry)
{
    const double tol_d = 0.000001;
    for (const auto& r_node : rGeometry) {
        if (std::abs(r_node.FastGetSolutionStepValue(DISTANCE)) < tol_d) {
            return true;
        }
    }
    return false;
}

bool EmbeddedMLSConstraintProcess::IsNegative(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        }
    }
    return (n_neg == rGeometry.PointsNumber());
}

MLSShapeFunctionsFunctionType EmbeddedMLSConstraintProcess::GetMLSShapeFunctionsFunction()
{
    switch (mpModelPart->GetProcessInfo()[DOMAIN_SIZE]) {
        case 2:
            switch (mMLSExtensionOperatorOrder) {
                case 1:
                    return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(rPoints, rX, h, rN);};
                case 2:
                    return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<2,2>(rPoints, rX, h, rN);};
                default:
                    KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
            }
        case 3:
            switch (mMLSExtensionOperatorOrder) {
                case 1:
                    return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(rPoints, rX, h, rN);};
                case 2:
                    return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<3,2>(rPoints, rX, h, rN);};
                default:
                    KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
            }
        default:
            KRATOS_ERROR << "Wrong domain size. MLS shape functions utility cannot be set.";
    }
}

void EmbeddedMLSConstraintProcess::SetNegativeNodeSupportCloud(
    const NodeType& rNegativeNode,
    PointerVector<NodeType>& rCloudNodes,
    PointerVector<NodeType>& rPositiveNeighborNodes,
    Matrix& rCloudCoordinates)
{
    // Find the positive side support cloud of nodes
    // Note that we use an unordered_set to ensure that these are unique
    NodesCloudSetType aux_set;
    std::vector<NodeType::Pointer> cur_layer_nodes;
    std::vector<NodeType::Pointer> prev_layer_nodes;
    const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

    // Find the negative nodes neighbors that are in the surrogate boundary
    // This is to find the first layer of positive nodes neighboring a negative one
    auto& no_const_node = const_cast<NodeType&>(rNegativeNode);
    auto& r_nod_neigh = no_const_node.GetValue(NEIGHBOUR_NODES);
    const std::size_t n_neigh = r_nod_neigh.size();
    for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
        auto& r_neigh = r_nod_neigh[i_neigh];
        if (r_neigh.Is(BOUNDARY)) {
            NodeType::Pointer p_node = &r_neigh;
            aux_set.insert(p_node);
            prev_layer_nodes.push_back(p_node);
        }
    }
    // Get vector of pointers to all positive nodes, the node is directly connected to
    const std::size_t n_positive_neigbors = aux_set.size();
    rPositiveNeighborNodes.resize(n_positive_neigbors);
    std::size_t i_neighbor = 0;
    for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
        rPositiveNeighborNodes(i_neighbor++) = *it_set;
    }

    // Add first layer neighbors to map
    // Note that we check the order of the MLS interpolation to add nodes from enough interior layers
    for (std::size_t i_layer = 0; i_layer < n_layers; ++i_layer) {
        for (auto& p_node : prev_layer_nodes) {
            auto& r_pos_node_neigh = p_node->GetValue(NEIGHBOUR_NODES);
            const std::size_t n_neigh = r_pos_node_neigh.size();
            for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
                auto& r_neigh = r_pos_node_neigh[i_neigh];
                if (r_neigh.Is(ACTIVE)) {
                    NodeType::Pointer p_neigh = &r_neigh;
                    p_neigh->Set(INTERFACE, true);
                    aux_set.insert(p_neigh);
                    cur_layer_nodes.push_back(p_neigh);
                }
            }
        }

        prev_layer_nodes = cur_layer_nodes;
        cur_layer_nodes.clear();
    }

    // Check that the current nodes are enough to perform the MLS calculation
    // If not sufficient, add the nodal neighbors of the current nodes to the cloud of points
    const std::size_t n_cloud_nodes_temp = aux_set.size();
    KRATOS_ERROR_IF(n_cloud_nodes_temp == 0) << "Degenerated case with no neighbours. Check if the node " << rNegativeNode.Id() << " belongs to an intersected and isolated element." << std::endl;
    if (n_cloud_nodes_temp < GetRequiredNumberOfPoints()) {
        NodesCloudSetType aux_extra_set(aux_set);
        for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
            auto& r_it_set_neighs = (*it_set)->GetValue(NEIGHBOUR_NODES);
            const std::size_t n_neigh = r_it_set_neighs.size();
            for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh) {
                auto& r_neigh = r_it_set_neighs[i_neigh];
                if (r_neigh.Is(ACTIVE)) {
                    NodeType::Pointer p_neigh = &r_neigh;
                    aux_extra_set.insert(p_neigh);
                }
            }
        }
        aux_set = aux_extra_set;
    }

    // Sort the obtained nodes by id
    const std::size_t n_cloud_nodes = aux_set.size();
    rCloudNodes.resize(n_cloud_nodes);
    std::size_t aux_i = 0;
    for (auto it_set = aux_set.begin(); it_set != aux_set.end(); ++it_set) {
        rCloudNodes(aux_i++) = *it_set;
    }
    std::sort(rCloudNodes.ptr_begin(), rCloudNodes.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});

    // Fill the coordinates matrix
    rCloudCoordinates.resize(n_cloud_nodes, 3);
    IndexPartition<std::size_t>(n_cloud_nodes).for_each(array_1d<double,3>(), [&rCloudNodes, &rCloudCoordinates](std::size_t iNode, array_1d<double,3>& rAuxCoordTLS){
        noalias(rAuxCoordTLS) = rCloudNodes[iNode].Coordinates();
        rCloudCoordinates(iNode, 0) = rAuxCoordTLS[0];
        rCloudCoordinates(iNode, 1) = rAuxCoordTLS[1];
        rCloudCoordinates(iNode, 2) = rAuxCoordTLS[2];
    });
}

void EmbeddedMLSConstraintProcess::SetNegativeNodeIntersectionPoints(
    const NodeType& rNegativeNode,
    PointerVector<NodeType>& rPositiveNeighborNodes,
    Matrix& rIntersectionPointsCoordinates)
{
    // Resize coordinate matrix
    const std::size_t n_positive_neighbors = rPositiveNeighborNodes.size();
    rIntersectionPointsCoordinates.resize(n_positive_neighbors, 3);

    // Use distance values of negative node and each directly neighboring positive node to calculate intersection point between the negative node and that positive node
    const double dist_neg_node = rNegativeNode.FastGetSolutionStepValue(DISTANCE);
    for (std::size_t i_node = 0; i_node < n_positive_neighbors; ++i_node) {
        array_1d<double,3> int_pt_coord;
        // Distance of positive node
        auto p_pos_node = rPositiveNeighborNodes(i_node);
        const double dist_pos_node = p_pos_node->FastGetSolutionStepValue(DISTANCE);
        // Calculate coordinates of intersection point using the distance values and coordinates of the negative node and the positive node
        const double int_pt_rel_location = std::abs(dist_neg_node/(dist_pos_node-dist_neg_node));
        int_pt_coord = rNegativeNode + int_pt_rel_location * (*p_pos_node - rNegativeNode);
        // Add intersection point to coordinate matrix of intersection points for the negative node
        rIntersectionPointsCoordinates(i_node, 0) = int_pt_coord[0];
        rIntersectionPointsCoordinates(i_node, 1) = int_pt_coord[1];
        rIntersectionPointsCoordinates(i_node, 2) = int_pt_coord[2];
    }
}

double EmbeddedMLSConstraintProcess::CalculateKernelRadius(
    const Matrix& rCloudCoordinates,
    const array_1d<double,3>& rOrigin)
{
    const std::size_t n_nodes = rCloudCoordinates.size1();
    const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
        return std::pow(rCloudCoordinates(I,0) - rOrigin(0),2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1),2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2),2);
    });
    return std::sqrt(squared_rad);
}

std::size_t EmbeddedMLSConstraintProcess::GetRequiredNumberOfPoints()
{
    const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
    switch (n_dim) {
        case 2:
            switch (mMLSExtensionOperatorOrder) {
                case 1:
                    return 3;
                case 2:
                    return 6;
                default:
                    KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
            }
        case 3:
            switch (mMLSExtensionOperatorOrder) {
                case 1:
                    return 4;
                case 2:
                    return 10;
                default:
                    KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
            }
        default:
            KRATOS_ERROR << "Wrong domain size.";
    }
}

}; // namespace Kratos.
