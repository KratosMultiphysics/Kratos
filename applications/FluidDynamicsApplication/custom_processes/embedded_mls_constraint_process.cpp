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
#include <unordered_set>

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "embedded_mls_constraint_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

/* Public functions *******************************************************/

EmbeddedMLSConstraintProcess::EmbeddedMLSConstraintProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
{
    // Validate input settings with defaults
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Retrieve the required model parts
    const std::string model_part_name = ThisParameters["model_part_name"].GetString();
    mpModelPart = &rModel.GetModelPart(model_part_name).GetSubModelPart("fluid_computational_model_part");

    // Set to which elements constraints are applied
    mApplyToAllNegativeCutNodes = ThisParameters["apply_to_all_negative_cut_nodes"].GetBool();

    // Set the order of the MLS extension operator used in the MLS shape functions utility
    mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

    // Set whether nodal distances will be modified to avoid levelset zeros
    mAvoidZeroDistances = ThisParameters["avoid_zero_distances"].GetBool();

    // Set which elements will not be active
    mDeactivateNegativeElements = ThisParameters["deactivate_negative_elements"].GetBool();
    mDeactivateIntersectedElements = ThisParameters["deactivate_intersected_elements"].GetBool();
}

void EmbeddedMLSConstraintProcess::Execute()
{
    // Declare extension operator map for negative nodes of split elements
    NodesCloudMapType ext_op_map;

    if (mAvoidZeroDistances) { ModifyDistances(); }

    // Set the required interface flags
    // NOTE: this will deactivate all negative as well as intersected elements
    SetInterfaceFlags();

    // Calculate the conforming extension basis
    CalculateConformingExtensionBasis(ext_op_map);

    // Reactivate intersected and negative elements as well as their nodes depending on the process settings
    ReactivateElementsAndNodes();

    // Apply extension constraints to negative nodes of split elements
    ApplyExtensionConstraints(ext_op_map);
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

void EmbeddedMLSConstraintProcess::CalculateConformingExtensionBasis(NodesCloudMapType& rExtensionOperatorMap)
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
                    const std::size_t found = rExtensionOperatorMap.count(&r_node);
                    if (r_node.IsNot(ACTIVE) && !found) {

                        // Get the current negative node neighbours cloud
                        Matrix cloud_nodes_coordinates;
                        PointerVector<NodeType> cloud_nodes;
                        SetNegativeNodeSupportCloud(r_node, cloud_nodes, cloud_nodes_coordinates);

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

                        auto ext_op_key_data = std::make_pair(&r_node, cloud_data_vector);
                        rExtensionOperatorMap.insert(ext_op_key_data);

                        KRATOS_WATCH("Added mls cloud for a negative node of a cut element");
                    }
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::ApplyExtensionConstraints(NodesCloudMapType& rExtensionOperatorMap)
{
    // Initialize counter of master slave constraints
    ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;

    // Define variables to constrain
    std::vector<const Variable<double>*> variables;
    const std::size_t n_dim = mpModelPart->GetProcessInfo()[DOMAIN_SIZE];
    switch (n_dim) {
        case 2:
            variables = {&VELOCITY_X, &VELOCITY_Y, &PRESSURE};
            break;
        case 3:
            variables = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z, &PRESSURE};
            break;
        default:
            KRATOS_ERROR << "Wrong domain size.";
    }

    // Loop through all negative nodes of split elements (slave nodes)
    for (const auto& r_slave : rExtensionOperatorMap) {
        const auto p_slave_node = r_slave.first;
        auto& r_ext_op_data = r_slave.second;

        // EVEN FASTER??
        DofPointerVectorType slave_dofs;
        DofPointerVectorType master_dofs;
        MatrixType dof_relation_matrix;
        VectorType dof_offset;

        const std::size_t n_var = n_dim + 1;
        dof_offset.resize(n_var, false);
        dof_relation_matrix.resize(n_var, r_ext_op_data.size()*n_var, false);

        std::size_t it_var = 0;
        for (const auto var : variables) {
            slave_dofs.push_back(p_slave_node->pGetDof(*var));
            dof_offset(it_var++) = 0.0;
        }

        std::size_t it_cloud = 0;
        for (auto& r_cloud_data : r_ext_op_data) {
            const auto p_cloud_node = r_cloud_data.first;
            const double cloud_node_N = r_cloud_data.second;

            std::size_t it_var = 0;
            for (const auto var : variables) {
                master_dofs.push_back(p_cloud_node->pGetDof(*var));

                for (std::size_t j = 0; j < n_var; j++) {
                    dof_relation_matrix(j, it_cloud*n_var+it_var) = (j == it_var) ? cloud_node_N : 0.0;
                }
                it_var++;
            }
            it_cloud++;
        }

        mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
        master_dofs, slave_dofs, dof_relation_matrix, dof_offset);


        //p_slave_node->Set(SLAVE);

        //WORKING VERSION - FAST!!
        /*for (const auto var : variables) {

            DofPointerVectorType slave_dofs;
            DofPointerVectorType master_dofs;
            MatrixType dof_relation_matrix;
            VectorType dof_offset;

            master_dofs.reserve(1);
            slave_dofs.push_back(p_slave_node->pGetDof(*var));
            master_dofs.reserve(r_ext_op_data.size());
            dof_offset.resize(1, false);
            dof_offset(0) = 0.0;
            dof_relation_matrix.resize(1, r_ext_op_data.size(), false);

            std::size_t it_cloud = 0;
            for (auto& r_cloud_data : r_ext_op_data) {
                const auto p_cloud_node = r_cloud_data.first;
                const double cloud_node_N = r_cloud_data.second;

                master_dofs.push_back(p_cloud_node->pGetDof(*var));
                dof_relation_matrix(0, it_cloud++) = cloud_node_N;
            }

            mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
            master_dofs, slave_dofs, dof_relation_matrix, dof_offset);
            KRATOS_WATCH(master_dofs.size());
        }*/
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

void EmbeddedMLSConstraintProcess::ReactivateElementsAndNodes()
{
    // Make intersected elements and their nodes ACTIVE
    if ( !mDeactivateIntersectedElements ) {
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
    }
    // Make negative elements and their nodes ACTIVE
    if ( !mDeactivateNegativeElements ) {
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

        const std::array<std::string,4> components = {"PRESSURE","VELOCITY_X","VELOCITY_Y","VELOCITY_Z"};

        // Set to zero and fix the DOFs in the remaining inactive nodes
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
            auto it_node = rNodes.begin() + i_node;
            if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
                for (std::size_t i_var = 0; i_var < components.size(); i_var++){
                    const auto& r_double_var = KratosComponents<Variable<double>>::Get(components[i_var]);
                    // Fix the nodal DOFs
                    it_node->Fix(r_double_var);
                    // Set to zero the nodal DOFs
                    it_node->FastGetSolutionStepValue(r_double_var) = 0.0;
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::ModifyDistances()
{
    auto& r_nodes = mpModelPart->Nodes();
    double tol_d = 1.0e-11;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
        auto it_node = r_nodes.begin() + i;
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
    Matrix& rCloudCoordinates)
{
    // Find the positive side support cloud of nodes
    // Note that we use an unordered_set to ensure that these are unique
    NodesCloudSetType aux_set;
    std::vector<NodeType::Pointer> cur_layer_nodes;
    std::vector<NodeType::Pointer> prev_layer_nodes;
    const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

    // Find the negative nodes neighbours that are in the surrogate boundary
    // This is to find the first layer of positive nodes neighbouring a negative one
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

    // Add first layer neighbours to map
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
    // If not sufficient, add the nodal neighbours of the current nodes to the cloud of points
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

double EmbeddedMLSConstraintProcess::CalculateKernelRadius(
    const Matrix& rCloudCoordinates,
    const array_1d<double,3>& rOrigin)
{
    const std::size_t n_nodes = rCloudCoordinates.size1();
    const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
        return std::pow(rCloudCoordinates(I,0) - rOrigin(0), 2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1), 2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2), 2);
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
