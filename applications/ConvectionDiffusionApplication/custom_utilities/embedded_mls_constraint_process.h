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

#if !defined(KRATOS_EMBEDDED_CONSTRAINT_UTILITIES )
#define  KRATOS_EMBEDDED_CONSTRAINT_UTILITIES

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "containers/pointer_vector.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "processes/process.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"

// Application includes

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

/// Short class definition.
class EmbeddedMLSConstraintProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    using IndexType = ModelPart::IndexType;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

    using MLSShapeFunctionsAndGradientsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)>;

    using ElementSizeFunctionType = std::function<double(const GeometryType&)>;

    using NodesCloudSetType = std::unordered_set<NodeType::Pointer, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    using CloudDataVectorType = DenseVector<std::pair<NodeType::Pointer, double>>;

    using NodesCloudMapType = std::unordered_map<NodeType::Pointer, CloudDataVectorType, SharedPointerHasher<NodeType::Pointer>, SharedPointerComparator<NodeType::Pointer>>;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of EmbeddedMLSConstraintProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedMLSConstraintProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    EmbeddedMLSConstraintProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
    {
        // Validate input settings with defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Retrieve the required model parts
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

        // Set whether nodal distances will be modified to avoid levelset zeros
        mAvoidZeroDistances = ThisParameters["avoid_zero_distances"].GetBool();

        // Set which elements will not be active
        mDeactivateNegativeElements = ThisParameters["deactivate_negative_elements"].GetBool();
        mDeactivateIntersectedElements = ThisParameters["deactivate_intersected_elements"].GetBool();
    };

    /// Destructor.
    virtual ~EmbeddedMLSConstraintProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
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

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "mls_extension_operator_order" : 1,
            "avoid_zero_distances" : true,
            "deactivate_negative_elements" : true,
            "deactivate_intersected_elements" : false
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "EmbeddedMLSConstraintProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmbeddedMLSConstraintProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;

    std::size_t mMLSExtensionOperatorOrder;

    bool mAvoidZeroDistances;

    bool mDeactivateNegativeElements;
    bool mDeactivateIntersectedElements;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateConformingExtensionBasis(NodesCloudMapType& rExtensionOperatorMap)
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
                    }
                }
            }
        }
    }

    void ApplyExtensionConstraints(NodesCloudMapType& rExtensionOperatorMap)
    {
        // Initialize counter of master slave constraints
        ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;
        // Get variable to constrain 
        const auto& r_var = KratosComponents<Variable<double>>::Get("TEMPERATURE");

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

    void SetInterfaceFlags()
    {
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

    void ReactivateElementsAndNodes()
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

    void ModifyDistances()
    {
        auto& r_nodes = mpModelPart->Nodes();
        double tol_d = 1.0e-12;

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

    bool IsSplit(const GeometryType& rGeometry)
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

    bool IsNegative(const GeometryType& rGeometry)
    {
        std::size_t n_neg = 0;
        for (const auto& r_node : rGeometry) {
            if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
                n_neg++;
            }
        }
        return (n_neg == rGeometry.PointsNumber());
    }

    MLSShapeFunctionsFunctionType GetMLSShapeFunctionsFunction()
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

    void SetNegativeNodeSupportCloud(
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

    double CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin)
    {
        const std::size_t n_nodes = rCloudCoordinates.size1();
        const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
            return std::pow(rCloudCoordinates(I,0) - rOrigin(0), 2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1), 2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2), 2);
        });
        return std::sqrt(squared_rad);
    }

    std::size_t GetRequiredNumberOfPoints()
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

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // EmbeddedMLSConstraintProcess& operator=(EmbeddedMLSConstraintProcess const& rOther);

    /// Copy constructor.
    //EmbeddedMLSConstraintProcess(EmbeddedMLSConstraintProcess const& rOther);

    ///@}

}; // Class EmbeddedMLSConstraintProcess

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_CONSTRAINT_UTILITIES  defined
