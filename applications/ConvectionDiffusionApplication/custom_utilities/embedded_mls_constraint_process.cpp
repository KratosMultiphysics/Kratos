// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
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
#include "utilities/reduction_utilities.h"

// Application includes
#include "embedded_mls_constraint_process.h"

namespace Kratos
{

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
    mpModelPart = &rModel.GetModelPart(model_part_name);

    // Retrieve the unknown variable
    mUnknownVariable = ThisParameters["unknown_variable"].GetString();

    // Set the order of the MLS extension operator used in the MLS shape functions utility
    mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

    // Set which elements will not be active
    mDeactivateNegativeElements = ThisParameters["deactivate_negative_elements"].GetBool();
    mDeactivateIntersectedElements = ThisParameters["deactivate_intersected_elements"].GetBool();
}

void EmbeddedMLSConstraintProcess::Execute()
{
    // Declare extension operator map for negative nodes of split elements
    NodesCloudMapType ext_op_map;

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
            // Find the intersected element negative nodes
            for (auto& r_node : r_geom) {
                // Check whether node is negative (not ACTIVE) or already was added because of a previous element
                const std::size_t found = rExtensionOperatorMap.count(&r_node);
                if (r_node.IsNot(ACTIVE) && !found) {

                    // Get the current negative node neighbors cloud
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
                }
            }
        }
    }
}

void EmbeddedMLSConstraintProcess::ApplyExtensionConstraints(NodesCloudMapType& rExtensionOperatorMap)
{
    // Initialize counter of master slave constraints
    ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;
    // Get variable to constrain
    const auto& r_var = KratosComponents<Variable<double>>::Get(mUnknownVariable);

    // Loop through all negative nodes of split elements (slave nodes)
    for (auto it_slave = rExtensionOperatorMap.begin(); it_slave != rExtensionOperatorMap.end(); ++it_slave) {
        auto p_slave_node = std::get<0>(*it_slave);
        auto& r_ext_op_data = rExtensionOperatorMap[p_slave_node];

        // Add one master slave constraint for every node of the support cloud (master) of the negative node (slave)
        // The contributions of each master will be summed up in the BuilderAndSolver to give an equation for the slave dof
        for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
            auto& r_node_data = *it_data;
            auto p_support_node = std::get<0>(r_node_data);
            const double support_node_N = std::get<1>(r_node_data);

            // Add master slave constraint, the support node MLS shape function value N serves as weight of the constraint
            mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
            *p_support_node, r_var, *p_slave_node, r_var,
            support_node_N, 0.0);
        }
    }
}

void EmbeddedMLSConstraintProcess::SetInterfaceFlags()
{
    // Initialize flags to false
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.Set(ACTIVE, false);        // Nodes that belong to completely positive elements
        rNode.Set(INTERFACE, false);     // Nodes that belong to the clouds points for MLS shape functions
    });
    block_for_each(mpModelPart->Elements(), [](Element& rElement){
        rElement.Set(ACTIVE, false);     // Elements in the positive distance region (the ones to be assembled)
        rElement.Set(INTERFACE, false);  // Positive distance elements owning the cloud nodes
    });

    // Find active regions
    for (auto& rElement : mpModelPart->Elements()) {
        auto& r_geom = rElement.GetGeometry();
        if (!IsSplit(r_geom) && !IsNegative(r_geom)) {
            // Mark the positive element as ACTIVE to assemble it
            rElement.Set(ACTIVE, true);
            // Mark the active element nodes as ACTIVE
            for (auto& rNode : r_geom) {
                rNode.Set(ACTIVE, true);
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

EmbeddedMLSConstraintProcess::MLSShapeFunctionsFunctionType EmbeddedMLSConstraintProcess::GetMLSShapeFunctionsFunction()
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
    NodesCloudSetType cloud_nodes_set;
    std::vector<NodeType::Pointer> current_layer_nodes;
    std::vector<NodeType::Pointer> previous_layer_nodes;
    const std::size_t n_layers = mMLSExtensionOperatorOrder + 1;

    // Add negative node to previous layer, in order to find the negative nodes neighbors that are in the surrogate boundary
    // This is to find the first layer of positive nodes neighboring a negative one
    NodeType* p_negative_node = &const_cast<NodeType&>(rNegativeNode);
    previous_layer_nodes.push_back(p_negative_node);

    // Note that we check the order of the MLS interpolation to add nodes from enough interior layers
    for (std::size_t i_layer = 0; i_layer < n_layers+1; ++i_layer) {
        for (auto& p_node : previous_layer_nodes) {
            auto& r_neighbors = p_node->GetValue(NEIGHBOUR_NODES);
            const std::size_t n_neighbors = r_neighbors.size();
            for (std::size_t i_neighbor = 0; i_neighbor < n_neighbors; ++i_neighbor) {
                auto& r_neighbor = r_neighbors[i_neighbor];
                if (r_neighbor.Is(ACTIVE)) {
                    NodeType::Pointer p_neighbor = &r_neighbor;
                    p_neighbor->Set(INTERFACE, true);
                    cloud_nodes_set.insert(p_neighbor);
                    current_layer_nodes.push_back(p_neighbor);
                }
            }
        }
        previous_layer_nodes = current_layer_nodes;
        current_layer_nodes.clear();
    }

    // Check that the current number of support nodes is enough to perform the MLS calculation
    // If not sufficient, add the nodal neighbors of the current nodes to the support cloud
    const std::size_t n_cloud_nodes_temp = cloud_nodes_set.size();
    KRATOS_ERROR_IF(n_cloud_nodes_temp == 0) << "Degenerated case with no neighbors. Check if the node " << rNegativeNode.Id() << " belongs to an intersected and isolated element." << std::endl;
    if (n_cloud_nodes_temp < GetRequiredNumberOfPoints()) {
        NodesCloudSetType cloud_nodes_set_extra(cloud_nodes_set);
        for (auto it_set = cloud_nodes_set.begin(); it_set != cloud_nodes_set.end(); ++it_set) {
            auto& r_neighbors = (*it_set)->GetValue(NEIGHBOUR_NODES);
            const std::size_t n_neighbors = r_neighbors.size();
            for (std::size_t i_neighbor = 0; i_neighbor < n_neighbors; ++i_neighbor) {
                auto& r_neighbor = r_neighbors[i_neighbor];
                if (r_neighbor.Is(ACTIVE)) {
                    NodeType::Pointer p_neighbor = &r_neighbor;
                    p_neighbor->Set(INTERFACE, true);
                    cloud_nodes_set_extra.insert(p_neighbor);
                }
            }
        }
        cloud_nodes_set = cloud_nodes_set_extra;
    }

    // Sort the obtained nodes by id
    const std::size_t n_cloud_nodes = cloud_nodes_set.size();
    rCloudNodes.resize(n_cloud_nodes);
    std::size_t aux_i = 0;
    for (auto it_set = cloud_nodes_set.begin(); it_set != cloud_nodes_set.end(); ++it_set) {
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
