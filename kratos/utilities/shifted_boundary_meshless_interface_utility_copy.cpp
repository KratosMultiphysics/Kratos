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
//

// System includes

// External includes

// Project includes
#include "containers/pointer_vector.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/element_size_calculator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/shifted_boundary_meshless_interface_utility_copy.h"

namespace Kratos
{

namespace
{
    using GeometryType = ShiftedBoundaryMeshlessInterfaceUtilityCopy::GeometryType;
    using ModifiedShapeFunctionsFactoryType = ShiftedBoundaryMeshlessInterfaceUtilityCopy::ModifiedShapeFunctionsFactoryType;

    /**
     * @brief Check if current geometry is split
     * This method checks if current geometry is split from the nodal historical level set values
     * @param rGeometry Geometry to be checked
     * @param rLevelSetVariable Variable storing the level set function
     * @return true If split
     * @return false If not split
     */
    bool IsSplit(
        const GeometryType& rGeometry,
        const Variable<double>& rLevelSetVariable)
    {
        std::size_t n_neg = 0;
        std::size_t n_pos = 0;
        for (const auto& r_node : rGeometry) {
            if (r_node.FastGetSolutionStepValue(rLevelSetVariable) < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }
        return (n_pos != 0 && n_neg != 0);
    }

    /**
     * @brief Check if current geometry is in negative side
     * This method checks if current geometry has negative level set value in all the nodes
     * @param rGeometry Geometry to be checked
     * @param rLevelSetVariable Variable storing the level set function
     * @return true If all nodes are negative
     * @return false If there is one or more positive nodes
     */
    bool IsNegative(
        const GeometryType& rGeometry,
        const Variable<double>& rLevelSetVariable)
    {
        std::size_t n_neg = 0;
        for (const auto& r_node : rGeometry) {
            if (r_node.FastGetSolutionStepValue(rLevelSetVariable) < 0.0) {
                n_neg++;
            }
        }
        return (n_neg == rGeometry.PointsNumber());
    }

    /**
     * @brief Set the nodal distances vector
     * This method saves the nodal historical values of level set in the provided vector
     * @param rGeometry Geometry from which the nodal values are retrieved
     * @param rLevelSetVariable Variable storing the level set function
     * @param rNodalDistances Vector container to store the distance values
     */
    void SetNodalDistancesVector(
        const GeometryType& rGeometry,
        const Variable<double>& rLevelSetVariable,
        Vector& rNodalDistances)
    {
        const std::size_t n_nodes = rGeometry.PointsNumber();
        if (rNodalDistances.size() != n_nodes) {
            rNodalDistances.resize(n_nodes);
        }
        std::size_t i = 0;
        for (const auto& r_node : rGeometry) {
            rNodalDistances[i++] = r_node.FastGetSolutionStepValue(rLevelSetVariable);
        }
    }

    /**
     * @brief Get the standard modified shape functions factory object
     * This function returns a prototype for the split shape functions calculation from the provided geometry
     * @param rGeometry Input geometry
     * @return ModifiedShapeFunctionsFactoryType Factory to be used for the split shape functions calculation
     */
    ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType::Pointer pGeometry, const Vector& rNodalDistances)->ModifiedShapeFunctions::UniquePointer{
                    return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }
}

    ShiftedBoundaryMeshlessInterfaceUtilityCopy::ShiftedBoundaryMeshlessInterfaceUtilityCopy(
        Model& rModel,
        Parameters ThisParameters)
    {
        // Validate input settings with defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Retrieve the required model parts
        const std::string model_part_name = ThisParameters["model_part_name"].GetString();
        const std::string boundary_sub_model_part_name = ThisParameters["boundary_sub_model_part_name"].GetString();
        mpModelPart = &rModel.GetModelPart(model_part_name);
        if (mpModelPart->HasSubModelPart(boundary_sub_model_part_name)) {
            mpBoundarySubModelPart = &(mpModelPart->GetSubModelPart(boundary_sub_model_part_name));
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfNodes() != 0) << "Provided SBM model part has nodes." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfElements() != 0) << "Provided SBM model part has elements." << std::endl;
            KRATOS_WARNING_IF("ShiftedBoundaryMeshlessInferfaceProcess", mpBoundarySubModelPart->NumberOfConditions() != 0) << "Provided SBM model part has conditions." << std::endl;
        } else {
            mpBoundarySubModelPart = &(mpModelPart->CreateSubModelPart(boundary_sub_model_part_name));
        }

        // Save a pointer to variable storing the level set
        const std::string levelset_variable_name = ThisParameters["levelset_variable_name"].GetString();
        if (KratosComponents<Variable<double>>::Has(levelset_variable_name)) {
            mpLevelSetVariable = &KratosComponents<Variable<double>>::Get(levelset_variable_name);
        } else {
            KRATOS_ERROR << "Provided 'levelset_variable_name' " << levelset_variable_name << " is not registered." << std::endl;
        }

        // Set the order of the MLS extension operator used in the MLS shape functions utility
        mMLSExtensionOperatorOrder = ThisParameters["mls_extension_operator_order"].GetInt();

        // If true, the basis is created such that it is conforming with the linear FE space of the surrogate boundary
        mConformingBasis = ThisParameters["conforming_basis"].GetBool();

        // If true, the basis is created such that the surrogate boundary gradient is kept
        const std::string ext_op_type = ThisParameters["extension_operator_type"].GetString();
        if (ext_op_type == "MLS") {
            mExtensionOperator = ExtensionOperator::MLS;
        } else if (ext_op_type == "RBF") {
            mExtensionOperator = ExtensionOperator::RBF;
        } else if (ext_op_type == "gradient_based") {
            mExtensionOperator = ExtensionOperator::GradientBased;
        } else {
            KRATOS_ERROR << "Wrong 'extension_operator_type' provided. Available options are 'MLS', 'RBF' and 'gradient_based'." << std::endl;
        }

        // Check that the basis settings are correct
        KRATOS_ERROR_IF(mExtensionOperator == ExtensionOperator::RBF && !mConformingBasis)
            << "'RBF' extension operator can only be used with conforming basis." << std::endl;
        KRATOS_ERROR_IF(mExtensionOperator == ExtensionOperator::GradientBased && !mConformingBasis)
            << "'gradient_based' extension operator can only be used with conforming basis." << std::endl;

        // Set the SBD contion prototype to be used in the condition creation
        std::string interface_condition_name = ThisParameters["sbm_interface_condition_name"].GetString();
        KRATOS_ERROR_IF(interface_condition_name == "") << "SBM interface condition has not been provided." << std::endl;
        mpConditionPrototype = &KratosComponents<Condition>::Get(interface_condition_name);
    };

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::CalculateExtensionOperator()
    {
        if (mExtensionOperator == ExtensionOperator::GradientBased) {
            CalculateGradientBasedConformingExtensionBasis();
        } else {
            if (mConformingBasis) {
                CalculateMeshlessBasedConformingExtensionBasis();
            } else {
                CalculateMeshlessBasedNonConformingExtensionBasis();
            }
        }
    }

    const Parameters ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "boundary_sub_model_part_name" : "",
            "sbm_interface_condition_name" : "",
            "conforming_basis" : true,
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1,
            "levelset_variable_name" : "DISTANCE"
        })" );

        return default_parameters;
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::CalculateGradientBasedConformingExtensionBasis()
    {
        // Get DOMAIN_SIZE
        const auto& r_process_info = mpModelPart->GetProcessInfo();
        KRATOS_ERROR_IF_NOT(r_process_info.Has(DOMAIN_SIZE)) << "'DOMAIN_SIZE' not present in ProcessInfo container." << std::endl;
        std::size_t domain_size = r_process_info[DOMAIN_SIZE];

        // Set the required interface flags
        // SetInterfaceFlags();             // OOOOOOOOOOCCHIOOOOOOOOO

        // Create the nodal gradient projection weights
        // This means to build a gradient reconstruction for each node in the surrogate boundary from the positive neighbour values
        auto sur_bd_nodes_map = SetSurrogateBoundaryNodalGradientWeights();

        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;
        const auto& r_begin_geom = mpModelPart->ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(r_begin_geom);

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Loop the elements to create the negative nodes basis
        NodesCloudMapType ext_op_map;
        for (auto& rElement : mpModelPart->Elements()) {

            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();
            if (IsSplit(*p_geom, *mpLevelSetVariable)) {
                // Find the intersected element negative nodes
                for (auto& r_node : *p_geom) {

                    // Check if the current node weight have been already computed
                    const std::size_t found = ext_op_map.count(&r_node);
                    if (r_node.IsNot(ACTIVE) && !found) {

                        // Get the neighbouring nodes that belong to the surrogate boundary
                        std::vector<std::size_t> bd_nodes_ids;
                        auto& r_neigh_vect = r_node.GetValue(NEIGHBOUR_NODES);
                        for (auto& r_bd_neigh : r_neigh_vect) {
                            if (r_bd_neigh.Is(BOUNDARY)) {
                                bd_nodes_ids.push_back(r_bd_neigh.Id());
                            }
                        }

                        // Calculate the weight of each surrogate node from its distance
                        // This is a weighted average between all the surrogate boundary nodes neighbouring the current extension node
                        double inv_tot_dist = 0.0;
                        std::size_t aux_i = 0;
                        Vector w_sur_nodes(bd_nodes_ids.size());
                        for (std::size_t bd_node_id : bd_nodes_ids) {
                            const auto p_sur_bd_node = mpModelPart->pGetNode(bd_node_id);
                            const double aux_dist = 1.0 / r_node.Distance(*p_sur_bd_node);
                            inv_tot_dist += aux_dist;
                            w_sur_nodes(aux_i) = aux_dist;
                            ++aux_i;
                        }
                        w_sur_nodes /= inv_tot_dist;

                        // Set the cloud weights
                        // In here we calculate the final weights for the extension operator
                        // Note that this includes the gradient projection from the surrogate nodes
                        aux_i = 0;
                        std::map<std::size_t, double> cloud_data_map;
                        for (std::size_t sur_bd_id : bd_nodes_ids) {
                            const double w_sur_node = w_sur_nodes(aux_i);
                            const auto p_sur_bd_node = mpModelPart->pGetNode(sur_bd_id);
                            const auto& r_sur_bd_dn_dx_data = sur_bd_nodes_map[sur_bd_id];
                            array_1d<double,3> proj_vect = r_node.Coordinates() - p_sur_bd_node->Coordinates();

                            // Add projected gradients weights
                            for (auto& dn_dx_data : r_sur_bd_dn_dx_data) {
                                // Calculate gradient edge projection
                                const auto& aux_dn_dx = dn_dx_data.second;
                                double dn_dx_proj_w = 0.0;
                                for (std::size_t d = 0; d < domain_size; ++d) {
                                    dn_dx_proj_w += aux_dn_dx[d] * proj_vect[d];
                                }
                                dn_dx_proj_w *= w_sur_node;

                                // Save current value
                                const std::size_t aux_id = dn_dx_data.first;
                                auto it_found = cloud_data_map.find(aux_id);
                                if (it_found != cloud_data_map.end()) {
                                    auto& r_val = it_found->second;
                                    r_val += dn_dx_proj_w;
                                } else {
                                    cloud_data_map.insert(std::make_pair(aux_id, dn_dx_proj_w));
                                }
                            }

                            // Add the current surrogate node weight
                            // Note that the "negative" node value is the surrogate nodes average one plus de nodal gradients average contribution
                            auto it_found = cloud_data_map.find(sur_bd_id);
                            if (it_found != cloud_data_map.end()) {
                                auto& r_val = it_found->second;
                                r_val += w_sur_node;
                            } else {
                                cloud_data_map.insert(std::make_pair(sur_bd_id, w_sur_node));
                            }

                            ++aux_i;
                        }

                        // Save the extension operator nodal data
                        std::size_t n_cl_nod = cloud_data_map.size();
                        double w_tot = 0.0;
                        CloudDataVectorType cloud_data_vector(n_cl_nod);
                        std::size_t i_cl_nod = 0;
                        for (auto& it_data : cloud_data_map) {
                            auto p_cl_node = mpModelPart->pGetNode(it_data.first);
                            auto i_data = std::make_pair(p_cl_node, it_data.second);
                            w_tot += it_data.second;
                            cloud_data_vector(i_cl_nod) = i_data;
                            i_cl_nod++;
                        }
                        KRATOS_ERROR_IF(std::abs(w_tot - 1.0 ) > 1.0e-12) << "Non-unit total weight " << w_tot << " for ext. operator in node " << r_node.Id() << std::endl;

                        auto ext_op_key_data = std::make_pair(&r_node, cloud_data_vector);
                        ext_op_map.insert(ext_op_key_data);
                    }
                }
            }
        }

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();
            if (IsSplit(*p_geom, *mpLevelSetVariable)) {
                // Set up the distances vector
                const auto& r_geom = *p_geom;
                const std::size_t n_nodes = r_geom.PointsNumber();
                Vector nodal_distances(n_nodes);
                SetNodalDistancesVector(r_geom, *mpLevelSetVariable, nodal_distances);

                // Set the modified shape functions pointer and calculate the positive interface data
                auto p_mod_sh_func = p_mod_sh_func_factory(p_geom, nodal_distances);
                Vector pos_int_w;
                Matrix pos_int_N;
                std::vector<array_1d<double,3>> pos_int_n;
                typename ModifiedShapeFunctions::ShapeFunctionsGradientsType pos_int_DN_DX;
                //TODO: Add a method without the interface gradients
                p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(pos_int_N, pos_int_DN_DX, pos_int_w, GeometryData::IntegrationMethod::GI_GAUSS_2);
                p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(pos_int_n, GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Calculate parent element size for the SBM BC imposition
                const double h = p_element_size_func(*p_geom);

                // Create an auxiliary set with all the cloud nodes that affect the current element
                NodesCloudSetType cloud_nodes_set;
                for (auto& r_node : r_geom) {
                    NodeType::Pointer p_node = &r_node;
                    if (r_node.Is(ACTIVE)) {
                        cloud_nodes_set.insert(p_node);
                    } else {
                        auto& r_ext_op_data = ext_op_map[p_node];
                        for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                            auto& p_node = std::get<0>(*it_data);
                            cloud_nodes_set.insert(p_node);
                        }
                    }
                }

                // Save previous resuls in a pointer vector to be used in the creation of the condition
                // Note that the obtained cloud is sorted by id to properly get the extension operator data
                PointerVector<NodeType> cloud_nodes_vector;
                const std::size_t n_cloud_nodes = cloud_nodes_set.size();
                cloud_nodes_vector.resize(n_cloud_nodes);
                std::size_t aux_i = 0;
                for (auto it_set = cloud_nodes_set.begin(); it_set != cloud_nodes_set.end(); ++it_set) {
                    cloud_nodes_vector(aux_i++) = *it_set;
                }
                std::sort(cloud_nodes_vector.ptr_begin(), cloud_nodes_vector.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});

                // Iterate the interface Gauss pts.
                DenseVector<double> i_g_N;
                DenseMatrix<double> i_g_DN_DX;
                array_1d<double,3> i_g_coords;
                const std::size_t n_int_pts = pos_int_w.size();
                for (std::size_t i_g = 0; i_g < n_int_pts; ++i_g) {
                    // Calculate Gauss pt. coordinates
                    i_g_N = row(pos_int_N, i_g);
                    i_g_DN_DX = pos_int_DN_DX[i_g];
                    noalias(i_g_coords) = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(i_g_coords) += i_g_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Initialize the extension operator containers
                    const std::size_t n_cl_nodes = cloud_nodes_vector.size();
                    const std::size_t n_dim = r_geom.WorkingSpaceDimension();
                    Vector N_container = ZeroVector(n_cl_nodes);
                    Matrix DN_DX_container = ZeroMatrix(n_cl_nodes, n_dim);

                    // Loop the nodes that are involved in the current element
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        const auto& r_node = r_geom[i_node];
                        if (r_node.Is(ACTIVE)) {
                            // If it is ACTIVE (positive) side add the standard shape function contribution
                            // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                            for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                                auto& p_cl_node = cloud_nodes_vector(i_cl);
                                if (r_node.Id() == p_cl_node->Id()) {
                                    N_container(i_cl) += i_g_N(i_node);
                                    for (std::size_t d = 0; d < n_dim; ++d) {
                                        DN_DX_container(i_cl,d) += i_g_DN_DX(i_node,d);
                                    }
                                    break;
                                }
                            }
                        } else {
                            // Get the weight as the corresponding nodal shape function value
                            const double i_node_N = i_g_N(i_node);
                            const auto i_node_grad_N = row(i_g_DN_DX, i_node);

                            // If it is not ACTIVE (negative side) search for the extension operator data
                            auto p_node = r_geom(i_node);
                            auto& ext_op_data = ext_op_map[p_node];

                            // Loop the current negative node to get its extrapolation operator data and apply the weight to make the basis conformant
                            // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                            for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                                auto& r_node_data = *it_data;
                                std::size_t data_node_id = (std::get<0>(r_node_data))->Id();
                                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                                    auto& p_cl_node = cloud_nodes_vector(i_cl);
                                    if (p_cl_node->Id() == data_node_id) {
                                        const double i_cl_node_N = std::get<1>(r_node_data);
                                        N_container(i_cl) += i_node_N * i_cl_node_N;
                                        for (std::size_t d = 0; d < n_dim; ++d) {
                                            DN_DX_container(i_cl,d) += i_node_grad_N(d) * i_cl_node_N;
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    double total_weights = 0.0;
                    for (std::size_t i = 0; i < cloud_nodes_vector.size(); ++i) {
                        total_weights += N_container(i);
                    }
                    KRATOS_ERROR_IF(std::abs(total_weights - 1.0) > 1.0e-12) << "Total weights is non-unit in split elem. " << rElement.Id() << std::endl;

                    // Create a new condition with a geometry made up with the basis nodes
                    auto p_prop = rElement.pGetProperties();
                    auto p_cond = mpConditionPrototype->Create(++max_cond_id, cloud_nodes_vector, p_prop);
                    p_cond->Set(ACTIVE, true);
                    mpBoundarySubModelPart->AddCondition(p_cond);

                    // Store the SBM BC data in the condition database
                    p_cond->SetValue(ELEMENT_H, h);
                    const double n_norm = norm_2(pos_int_n[i_g]);
                    p_cond->SetValue(NORMAL, pos_int_n[i_g] / n_norm);
                    p_cond->SetValue(INTEGRATION_WEIGHT, pos_int_w[i_g]);
                    p_cond->SetValue(INTEGRATION_COORDINATES, i_g_coords);
                    p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
                    p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);
                }
            }
        }
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::CalculateMeshlessBasedConformingExtensionBasis()
    {
        // Set the required interface flags
        // SetInterfaceFlags();             // OOOOOOOOOOCCHIOOOOOOOOO

        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;
        const auto& r_begin_geom = mpModelPart->ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the MLS shape functions function
        auto p_meshless_sh_func =
        mExtensionOperator == ExtensionOperator::MLS ? GetMLSShapeFunctionsFunction() : GetRBFShapeFunctionsFunction();

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(r_begin_geom);

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Loop the elements to create the negative nodes extension basis
        NodesCloudMapType ext_op_map;
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();
            if (IsSplit(*p_geom, *mpLevelSetVariable)) {
                // Find the intersected element negative nodes
                for (auto& r_node : *p_geom) {
                    const std::size_t found = ext_op_map.count(&r_node);
                    if (r_node.IsNot(ACTIVE) && !found) {
                        // Get the current negative node neighbours cloud
                        Matrix cloud_nodes_coordinates;
                        PointerVector<NodeType> cloud_nodes;
                        SetNegativeNodeSupportCloud(r_node, cloud_nodes, cloud_nodes_coordinates);

                        // Calculate the extension basis in the current negative node
                        Vector N_container;
                        const array_1d<double,3> r_coords = r_node.Coordinates();
                        const double kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
                        p_meshless_sh_func(cloud_nodes_coordinates, r_coords, kernel_rad, N_container);

                        // Save the extension operator nodal data
                        std::size_t n_cl_nod = cloud_nodes.size();
                        CloudDataVectorType cloud_data_vector(n_cl_nod);
                        for (std::size_t i_cl_nod = 0; i_cl_nod < n_cl_nod; ++i_cl_nod) {
                            auto p_cl_node = cloud_nodes(i_cl_nod);
                            auto i_data = std::make_pair(p_cl_node, N_container[i_cl_nod]);
                            cloud_data_vector(i_cl_nod) = i_data;
                        }

                        auto ext_op_key_data = std::make_pair(&r_node, cloud_data_vector);
                        ext_op_map.insert(ext_op_key_data);
                    }
                }
            }
        }

        // Create the interface conditions
        //TODO: THIS CAN BE PARALLEL (WE JUST NEED TO MAKE CRITICAL THE CONDITION ID UPDATE)
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();
            if (IsSplit(*p_geom, *mpLevelSetVariable)) {
                // Set up the distances vector
                const auto& r_geom = *p_geom;
                const std::size_t n_nodes = r_geom.PointsNumber();
                Vector nodal_distances(n_nodes);
                SetNodalDistancesVector(r_geom, *mpLevelSetVariable, nodal_distances);

                // Set the modified shape functions pointer and calculate the positive interface data
                auto p_mod_sh_func = p_mod_sh_func_factory(p_geom, nodal_distances);
                Vector pos_int_w;
                Matrix pos_int_N;
                std::vector<array_1d<double,3>> pos_int_n;
                typename ModifiedShapeFunctions::ShapeFunctionsGradientsType pos_int_DN_DX;
                //TODO: Add a method without the interface gradients
                p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(pos_int_N, pos_int_DN_DX, pos_int_w, GeometryData::IntegrationMethod::GI_GAUSS_2);
                p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(pos_int_n, GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Calculate parent element size for the SBM BC imposition
                const double h = p_element_size_func(*p_geom);

                // Create an auxiliary set with all the cloud nodes that affect the current element
                NodesCloudSetType cloud_nodes_set;
                for (auto& r_node : r_geom) {
                    NodeType::Pointer p_node = &r_node;
                    if (r_node.Is(ACTIVE)) {
                        cloud_nodes_set.insert(p_node);
                    } else {
                        auto& r_ext_op_data = ext_op_map[p_node];
                        for (auto it_data = r_ext_op_data.begin(); it_data != r_ext_op_data.end(); ++it_data) {
                            auto& p_node = std::get<0>(*it_data);
                            cloud_nodes_set.insert(p_node);
                        }
                    }
                }

                // Save previous resuls in a pointer vector to be used in the creation of the condition
                // Note that the obtained cloud is sorted by id to properly get the extension operator data
                PointerVector<NodeType> cloud_nodes_vector;
                const std::size_t n_cloud_nodes = cloud_nodes_set.size();
                cloud_nodes_vector.resize(n_cloud_nodes);
                std::size_t aux_i = 0;
                for (auto it_set = cloud_nodes_set.begin(); it_set != cloud_nodes_set.end(); ++it_set) {
                    cloud_nodes_vector(aux_i++) = *it_set;
                }
                std::sort(cloud_nodes_vector.ptr_begin(), cloud_nodes_vector.ptr_end(), [](NodeType::Pointer& pNode1, NodeType::Pointer rNode2){return (pNode1->Id() < rNode2->Id());});

                // Iterate the interface Gauss pts.
                DenseVector<double> i_g_N;
                DenseMatrix<double> i_g_DN_DX;
                array_1d<double,3> i_g_coords;
                const std::size_t n_int_pts = pos_int_w.size();
                for (std::size_t i_g = 0; i_g < n_int_pts; ++i_g) {
                    // Calculate Gauss pt. coordinates
                    i_g_N = row(pos_int_N, i_g);
                    i_g_DN_DX = pos_int_DN_DX[i_g];
                    noalias(i_g_coords) = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(i_g_coords) += i_g_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Initialize the extension operator containers
                    const std::size_t n_cl_nodes = cloud_nodes_vector.size();
                    const std::size_t n_dim = r_geom.WorkingSpaceDimension();
                    Vector N_container = ZeroVector(n_cl_nodes);
                    Matrix DN_DX_container = ZeroMatrix(n_cl_nodes, n_dim);

                    // Loop the nodes that are involved in the current element
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        const auto& r_node = r_geom[i_node];
                        if (r_node.Is(ACTIVE)) {
                            // If it is ACTIVE (positive) side add the standard shape function contribution
                            // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                            for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                                auto& p_cl_node = cloud_nodes_vector(i_cl);
                                if (r_node.Id() == p_cl_node->Id()) {
                                    N_container(i_cl) += i_g_N(i_node);
                                    for (std::size_t d = 0; d < n_dim; ++d) {
                                        DN_DX_container(i_cl,d) += i_g_DN_DX(i_node,d);
                                    }
                                    break;
                                }
                            }
                        } else {
                            // Get the weight as the corresponding nodal shape function value
                            const double i_node_N = i_g_N(i_node);
                            const auto i_node_grad_N = row(i_g_DN_DX, i_node);

                            // If it is not ACTIVE (negative side) search for the extension operator data
                            auto p_node = r_geom(i_node);
                            auto& ext_op_data = ext_op_map[p_node];

                            // Loop the current negative node to get its extrapolation operator data and apply the weight to make the basis conformant
                            // Note that we need to check for the ids to match in the geometry as nodes in the map are mixed
                            for (auto it_data = ext_op_data.begin(); it_data != ext_op_data.end(); ++it_data) {
                                auto& r_node_data = *it_data;
                                std::size_t data_node_id = (std::get<0>(r_node_data))->Id();
                                for (std::size_t i_cl = 0; i_cl < n_cl_nodes; ++i_cl) {
                                    auto& p_cl_node = cloud_nodes_vector(i_cl);
                                    if (p_cl_node->Id() == data_node_id) {
                                        const double i_cl_node_N = std::get<1>(r_node_data);
                                        N_container(i_cl) += i_node_N * i_cl_node_N;
                                        for (std::size_t d = 0; d < n_dim; ++d) {
                                            DN_DX_container(i_cl,d) += i_node_grad_N(d) * i_cl_node_N;
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    // Create a new condition with a geometry made up with the basis nodes
                    auto p_prop = rElement.pGetProperties();
                    auto p_cond = mpConditionPrototype->Create(++max_cond_id, cloud_nodes_vector, p_prop);
                    p_cond->Set(ACTIVE, true);
                    mpBoundarySubModelPart->AddCondition(p_cond);

                    // Store the SBM BC data in the condition database
                    p_cond->SetValue(ELEMENT_H, h);
                    const double n_norm = norm_2(pos_int_n[i_g]);
                    p_cond->SetValue(NORMAL, pos_int_n[i_g] / n_norm);
                    p_cond->SetValue(INTEGRATION_WEIGHT, pos_int_w[i_g]);
                    p_cond->SetValue(INTEGRATION_COORDINATES, i_g_coords);
                    p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
                    p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);
                }
            }
        }
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::CalculateMeshlessBasedNonConformingExtensionBasis()
    {
        // Set the required interface flags
        // SetInterfaceFlags();             // OOOOOOOOOOCCHIOOOOOOOOO

        // Set the modified shape functions factory
        // Note that unique geometry in the mesh is assumed
        KRATOS_ERROR_IF_NOT(mpModelPart->NumberOfElements())
            << "There are no elements in background mesh model part '" << mpModelPart->FullName() << "'." << std::endl;
        const auto& r_begin_geom = mpModelPart->ElementsBegin()->GetGeometry();
        auto p_mod_sh_func_factory = GetStandardModifiedShapeFunctionsFactory(r_begin_geom);

        // Get the MLS shape functions function
        auto p_mls_sh_func = GetMLSShapeFunctionsAndGradientsFunction();

        // Get the element size calculation function
        // Note that unique geometry in the mesh is assumed
        auto p_element_size_func = GetElementSizeFunction(r_begin_geom);

        // Get max condition id
        std::size_t max_cond_id = block_for_each<MaxReduction<std::size_t>>(mpModelPart->Conditions(), [](const Condition& rCondition){return rCondition.Id();});

        // Loop the elements to find the intersected ones
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto p_geom = rElement.pGetGeometry();

            // If the element is split, get the interface Gauss pt. data
            // If the element is negative, it is deactivated to not be assembled
            if (IsSplit(*p_geom, *mpLevelSetVariable)) {
                // Set the meshless cloud of point support for the MLS
                Matrix cloud_nodes_coordinates;
                PointerVector<NodeType> cloud_nodes;
                SetSplitElementSupportCloud(rElement, cloud_nodes, cloud_nodes_coordinates);

                // Set up the distances vector
                const auto& r_geom = *p_geom;
                const std::size_t n_nodes = r_geom.PointsNumber();
                Vector nodal_distances(n_nodes);
                SetNodalDistancesVector(r_geom, *mpLevelSetVariable, nodal_distances);

                // Set the modified shape functions pointer and calculate the positive interface data
                auto p_mod_sh_func = p_mod_sh_func_factory(p_geom, nodal_distances);
                Vector pos_int_w;
                Matrix pos_int_N;
                std::vector<array_1d<double,3>> pos_int_n;
                typename ModifiedShapeFunctions::ShapeFunctionsGradientsType pos_int_DN_DX;
                //TODO: Add a method without the interface gradients
                p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(pos_int_N, pos_int_DN_DX, pos_int_w, GeometryData::IntegrationMethod::GI_GAUSS_2);
                p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(pos_int_n, GeometryData::IntegrationMethod::GI_GAUSS_2);

                // Calculate parent element size for the SBM BC imposition
                const double h = p_element_size_func(*p_geom);

                // Iterate the interface Gauss pts.
                DenseVector<double> i_g_N;
                array_1d<double,3> i_g_coords;
                const std::size_t n_int_pts = pos_int_w.size();
                for (std::size_t i_g = 0; i_g < n_int_pts; ++i_g) {
                    // Calculate Gauss pt. coordinates
                    i_g_N = row(pos_int_N, i_g);
                    noalias(i_g_coords) = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(i_g_coords) += i_g_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Create a new condition with a geometry made up with the basis nodes
                    auto p_prop = rElement.pGetProperties();
                    auto p_cond = mpConditionPrototype->Create(++max_cond_id, cloud_nodes, p_prop);
                    p_cond->Set(ACTIVE, true);
                    mpBoundarySubModelPart->AddCondition(p_cond);

                    // Store the SBM BC data in the condition database
                    p_cond->SetValue(ELEMENT_H, h);
                    const double n_norm = norm_2(pos_int_n[i_g]);
                    p_cond->SetValue(NORMAL, pos_int_n[i_g] / n_norm);
                    p_cond->SetValue(INTEGRATION_WEIGHT, pos_int_w[i_g]);
                    p_cond->SetValue(INTEGRATION_COORDINATES, i_g_coords);

                    // Calculate the MLS shape functions and gradients and save in the database
                    Vector N_container;
                    Matrix DN_DX_container;
                    const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, i_g_coords);
                    p_mls_sh_func(cloud_nodes_coordinates, i_g_coords, mls_kernel_rad, N_container, DN_DX_container);
                    p_cond->SetValue(SHAPE_FUNCTIONS_VECTOR, N_container);
                    p_cond->SetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX, DN_DX_container);
                }
            }
        }
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::SetInterfaceFlags()
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
            if (IsSplit(r_geom, *mpLevelSetVariable)) {
                // Mark the intersected elements as BOUNDARY
                // The intersected elements are also flagged as non ACTIVE to avoid assembling them
                // Note that the split elements BC is applied by means of the extension operators
                rElement.Set(ACTIVE, false);
                rElement.Set(BOUNDARY, true);
                // Mark the positive intersected nodes as BOUNDARY
                for (auto& rNode : r_geom) {
                    if (!(rNode.FastGetSolutionStepValue(*mpLevelSetVariable) < 0.0)) {
                        rNode.Set(BOUNDARY, true);
                    }
                }
            } else if (IsNegative(r_geom, *mpLevelSetVariable)) {
                // Mark the negative element as non ACTIVE to avoid assembling it
                rElement.Set(ACTIVE, false);
            } else {
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
                    // The nodes will be flagged if required (MLS basis) when creating the cloud
                    auto p_neigh_elem = r_neigh_elems(i_face).get();
                    if (p_neigh_elem != nullptr) {
                        if (p_neigh_elem->Is(ACTIVE)) {
                            p_neigh_elem->Set(INTERFACE, true);
                        }
                    }
                }
            }
        }
    }

    ShiftedBoundaryMeshlessInterfaceUtilityCopy::MLSShapeFunctionsAndGradientsFunctionType ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetMLSShapeFunctionsAndGradientsFunction() const
    {
        switch (mpModelPart->GetProcessInfo()[DOMAIN_SIZE]) {
            case 2:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(rPoints, rX, h, rN, rDN_DX);};
                    case 2:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,2>(rPoints, rX, h, rN, rDN_DX);};
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            case 3:
                switch (mMLSExtensionOperatorOrder) {
                    case 1:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,1>(rPoints, rX, h, rN, rDN_DX);};
                    case 2:
                        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, Matrix& rDN_DX){
                            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,2>(rPoints, rX, h, rN, rDN_DX);};
                    default:
                        KRATOS_ERROR << "Wrong MLS extension operator order. Only linear (1) and quadratic (2) are supported.";
                }
            default:
                KRATOS_ERROR << "Wrong domain size. MLS shape functions utility cannot be set.";
        }
    }

    ShiftedBoundaryMeshlessInterfaceUtilityCopy::MeshlessShapeFunctionsFunctionType ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetMLSShapeFunctionsFunction() const
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

    ShiftedBoundaryMeshlessInterfaceUtilityCopy::MeshlessShapeFunctionsFunctionType ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetRBFShapeFunctionsFunction() const
    {
        return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
            RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, h, rN);
        };
    }

    ShiftedBoundaryMeshlessInterfaceUtilityCopy::ElementSizeFunctionType ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetElementSizeFunction(const GeometryType& rGeometry)
    {
        switch (rGeometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);};
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return [](const GeometryType& rGeometry)->double{return ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);};
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::SetSplitElementSupportCloud(
        const Element& rSplitElement,
        PointerVector<NodeType>& rCloudNodes,
        Matrix& rCloudCoordinates)
    {
        // Find the positive side support cloud of nodes
        // Note that we use an unordered_set to ensure that these are unique
        const auto& r_geometry = rSplitElement.GetGeometry();
        const std::size_t n_nodes = r_geometry.PointsNumber();
        const std::size_t n_layers = mMLSExtensionOperatorOrder;
        NodesCloudSetType aux_set;
        std::vector<NodeType::Pointer> cur_layer_nodes;
        std::vector<NodeType::Pointer> prev_layer_nodes;
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_geometry(i_node);
            if (p_node->Is(BOUNDARY)) {
                // Add current positive node to map
                aux_set.insert(p_node);
                prev_layer_nodes.push_back(p_node);

                // Add positive neighbours to map
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

                prev_layer_nodes.clear();
            }
        }

        // Check that the current nodes are enough to perform the MLS calculation
        // If not sufficient, add the nodal neighbours of the current nodes to the cloud of points
        const std::size_t n_cloud_nodes_temp = aux_set.size();
        KRATOS_ERROR_IF(n_cloud_nodes_temp == 0) << "Degenerated case with no neighbours. Check if the element " << rSplitElement.Id() << " is intersected and isolated." << std::endl;
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
        array_1d<double,3> aux_coord;
        rCloudCoordinates.resize(n_cloud_nodes, 3);
        for (std::size_t i_node = 0; i_node < n_cloud_nodes; ++i_node) {
            noalias(aux_coord) = rCloudNodes[i_node].Coordinates();
            rCloudCoordinates(i_node, 0) = aux_coord[0];
            rCloudCoordinates(i_node, 1) = aux_coord[1];
            rCloudCoordinates(i_node, 2) = aux_coord[2];
        }
    }

    void ShiftedBoundaryMeshlessInterfaceUtilityCopy::SetNegativeNodeSupportCloud(
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

    double ShiftedBoundaryMeshlessInterfaceUtilityCopy::CalculateKernelRadius(
        const Matrix& rCloudCoordinates,
        const array_1d<double,3>& rOrigin)
    {
        const std::size_t n_nodes = rCloudCoordinates.size1();
        const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
            return std::pow(rCloudCoordinates(I,0) - rOrigin(0),2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1),2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2),2);
        });
        return std::sqrt(squared_rad);
    }

    std::size_t ShiftedBoundaryMeshlessInterfaceUtilityCopy::GetRequiredNumberOfPoints()
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
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    std::unordered_map<std::size_t, std::map<std::size_t, Vector>> ShiftedBoundaryMeshlessInterfaceUtilityCopy::SetSurrogateBoundaryNodalGradientWeights()
    {
        std::unordered_map<std::size_t, std::map<std::size_t, Vector>> sur_bd_nodes_map;
        for (auto& r_bd_node : mpModelPart->Nodes()) {
            auto it_found = sur_bd_nodes_map.find(r_bd_node.Id());
            if (it_found == sur_bd_nodes_map.end() && r_bd_node.Is(BOUNDARY)) {

                // Set an auxilary map to calculate the current node nodal gradient contributions
                std::map<std::size_t, Vector> neigh_dn_dx_map;

                // Calculate the nodal gradient coefficients in the ACTIVE elements neighbouring current node
                double w_total = 0.0;
                int count = 0 ;

                auto& r_elem_neigh_vect = r_bd_node.GetValue(NEIGHBOUR_ELEMENTS);
                // KRATOS_WATCH(r_elem_neigh_vect.size())
                if (r_elem_neigh_vect.size() > 3 || r_elem_neigh_vect.size() < 3) { // Standard case
                    for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                        auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                        if (p_elem_neigh->IsNot(INTERFACE)) {
                            if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                                count = count + 1 ;
                                // Calculate the current element weight
                                const auto& r_geom = p_elem_neigh->GetGeometry();
                                const std::size_t n_nodes = r_geom.PointsNumber();
                                const double w = r_geom.DomainSize() / n_nodes;
                                w_total += w;

                                // Calculate the current element nodal gradient
                                // Note that in here we assume that the gradient is constant within the element (simplex geometries)
                                ShapeFunctionsGradientsType grad_N;
                                r_geom.ShapeFunctionsIntegrationPointsGradients(grad_N, GeometryData::IntegrationMethod::GI_GAUSS_1);

                                // Save the current element nodal weights
                                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                                    std::size_t aux_id = r_geom[i_node].Id();
                                    const Vector r_node_dn_dx = row(grad_N[0], i_node);
                                    auto it_found = neigh_dn_dx_map.find(aux_id);
                                    if (it_found != neigh_dn_dx_map.end()) {
                                        auto& r_val = it_found->second;
                                        r_val += w*r_node_dn_dx;
                                    } else {
                                        neigh_dn_dx_map.insert(std::make_pair(aux_id, w*r_node_dn_dx));
                                    }
                                }

                            }
                        }
                    }
                } else { 
                    // //////////////////////////////////////////////////////////// Se non va cancellare DA qua
                    // // r_elem_neigh_vect.size() == 3 --> uno solo disponibile
                    // // C'e' solo un elemento disponibile --> Solo due nodi per calcolare il gradiente
                    // // Aggiungiamo elementi:
                    // // 1.1 ottenere i due nodi vicini
                    // auto& r_node_neigh_vect = r_bd_node.GetValue(NEIGHBOUR_NODES);
                    // int lasciapassare = 0 ;
                    // auto p_node_neigh1 = r_node_neigh_vect(1).get(); // Initialize pointers
                    // auto p_node_neigh2 = r_node_neigh_vect(1).get();  // Initialize pointers
                    // for (std::size_t i_neigh_node = 0; i_neigh_node < r_node_neigh_vect.size(); ++i_neigh_node) {
                    //     auto p_node_neigh = r_node_neigh_vect(i_neigh_node).get();
                    //     double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                    //     if (phi_distance > 0 && p_node_neigh-> IsNot(BOUNDARY)) {
                    //         if (lasciapassare == 0) {
                    //             p_node_neigh1 = p_node_neigh ; 
                    //             lasciapassare = 1 ;
                    //         }else{
                    //             p_node_neigh2 = p_node_neigh ;
                    //         }
                             
                    //     }
                    // } 
                    // // Passare da Pointers a nodi veri e propri
                    // auto& r_node1 = *p_node_neigh1 ; 
                    // auto& r_node2 = *p_node_neigh2 ;
                    // // Abbiamo i due nodi, Possiamo cercare gli elementi vicini 
                    // auto& r_elem_neigh_vect1 = r_node1.GetValue(NEIGHBOUR_ELEMENTS);
                    // auto& r_elem_neigh_vect2 = r_node2.GetValue(NEIGHBOUR_ELEMENTS);
                    // // Ora aggiungiamo gli elementi a vettore di elementi già esistente e usiamo Unique
                    // // push_back
                    // for(auto& it : r_elem_neigh_vect1.GetContainer()) {
                    //     if (it->IsNot(BOUNDARY)){
                    //         r_elem_neigh_vect.push_back(it);
                    //     }         
                    // }     
                    // for(auto& it : r_elem_neigh_vect2.GetContainer()) {
                    //     if (it->IsNot(BOUNDARY)){
                    //         r_elem_neigh_vect.push_back(it);
                    //     }             
                    // }
                    // // Unique
                    // r_elem_neigh_vect.Unique() ;    
                    // // shrink_to_fit
                    // r_elem_neigh_vect.shrink_to_fit() ;
                    //////////////////////////////////////////////////////////// Se non va cancellare FINO A qua
                    for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                        auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                        if (p_elem_neigh->IsNot(INTERFACE)) {
                            if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                                count = count + 1 ;
                                // Calculate the current element weight
                                const auto& r_geom = p_elem_neigh->GetGeometry();
                                const std::size_t n_nodes = r_geom.PointsNumber();
                                // const double w = r_geom.DomainSize() / n_nodes;
                                auto center_coord = r_geom.Center();
                                const double w = 1/( (center_coord.X()-r_bd_node.X())*(center_coord.X()-r_bd_node.X()) + (center_coord.Y()-r_bd_node.Y())*(center_coord.Y()-r_bd_node.Y()) ) ;
                                w_total += w;

                                // Calculate the current element nodal gradient
                                // Note that in here we assume that the gradient is constant within the element (simplex geometries)
                                ShapeFunctionsGradientsType grad_N;
                                r_geom.ShapeFunctionsIntegrationPointsGradients(grad_N, GeometryData::IntegrationMethod::GI_GAUSS_1);

                                // Save the current element nodal weights
                                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                                    std::size_t aux_id = r_geom[i_node].Id();
                                    const Vector r_node_dn_dx = row(grad_N[0], i_node);
                                    auto it_found = neigh_dn_dx_map.find(aux_id);
                                    if (it_found != neigh_dn_dx_map.end()) {
                                        auto& r_val = it_found->second;
                                        r_val += w*r_node_dn_dx;
                                    } else {
                                        neigh_dn_dx_map.insert(std::make_pair(aux_id, w*r_node_dn_dx));
                                    }
                                }

                            }
                        }
                    }
                }
                // KRATOS_WATCH(count)
                // if (count == 0) {
                //     // KRATOS_ERROR << "Eccone uno" ;
                //     // KRATOS_WATCH("Eccone uno senza gradiente")
                // }

                // Divide current coefficients by the total weight
                for (auto& neigh_data : neigh_dn_dx_map) {
                    auto& r_val = neigh_data.second;
                    r_val /= w_total;
                }

                // Update the surrogate nodes data
                sur_bd_nodes_map.insert(std::make_pair(r_bd_node.Id(), neigh_dn_dx_map));
            }
        }
        return sur_bd_nodes_map;
    }

    
    //________________________________________________________________________________________________________________
    //________________________________________________________________________________________________________________
    //__________________________________NodalGradientWeightsForLonelyNodes_______________________________________
    std::unordered_map<std::size_t, std::map<std::size_t, Vector>> ShiftedBoundaryMeshlessInterfaceUtilityCopy::NodalGradientWeightsForLonelyNodes()
    {
        std::unordered_map<std::size_t, std::map<std::size_t, Vector>> sur_bd_nodes_map;
        for (auto& r_bd_node : mpModelPart->Nodes()) {
            auto it_found = sur_bd_nodes_map.find(r_bd_node.Id());
            if (it_found == sur_bd_nodes_map.end() && r_bd_node.Is(BOUNDARY)) {
                // Calculate the nodal gradient coefficients in the MARKER elements neighbouring current node
                // There are available elements? Let's count them
                int count = 0 ;        
                // Get the neighbours elements 
                auto& r_elem_neigh_vect = r_bd_node.GetValue(NEIGHBOUR_ELEMENTS);
                // Count how many of them are available for MPC
                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                    if (p_elem_neigh->IsNot(INTERFACE)) {
                        if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                            count = count + 1 ; 
                        }
                    }
                }
                KRATOS_WATCH(count)
                KRATOS_WATCH(r_bd_node.Id())
                // The node is either PROBLEMATIC or VERY_PROBLEMATIC
                if (count == 0) { // There are no available elements around the surrogate node
                    // KRATOS_WATCH(r_bd_node.Id())
                    // KRATOS_WATCH("Eccone uno senza <gradiente> !! ")
                    // Get the neighbours nodes of this surrogate node
                    auto& r_node_neigh_vect = r_bd_node.GetValue(NEIGHBOUR_NODES);
                    // Need to add a if loop for check that the node is not a VERY_PROBLEMATIC one
                    int count_neigh_nodes = 0 ; // Count the neighbours nodes that have DISTANCE > 0 
                    for (std::size_t i_neigh = 0; i_neigh < r_node_neigh_vect.size(); ++i_neigh) {
                        auto p_node_neigh = r_node_neigh_vect(i_neigh).get();
                        double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                        if (phi_distance > 0) {
                            count_neigh_nodes = count_neigh_nodes + 1 ;
                        }
                    }
                    KRATOS_WATCH(count_neigh_nodes)
                    if (count_neigh_nodes < 3) {
                        KRATOS_WATCH("VERY_PROBLEMATIC node")
                        if (count_neigh_nodes < 2) {
                            KRATOS_WATCH("DAMN PROBLEMATIC node ?")
                            exit(0) ;
                        }
                        int temp = 0 ; // Variable to find the two admissible neigh nodes
                        auto p_node_neigh1 = r_node_neigh_vect(1).get();
                        auto p_node_neigh2 = r_node_neigh_vect(1).get(); 
                        if (r_node_neigh_vect.size() == 0 ) {
                            KRATOS_WATCH("r_node_neigh_vect.size() == 0 !!! STOP") 
                            exit(0) ;
                        }
                        for (std::size_t i_neigh = 0; i_neigh < r_node_neigh_vect.size(); ++i_neigh) {
                            auto p_node_neigh = r_node_neigh_vect(i_neigh).get();
                            double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                            if (phi_distance > 0) { 
                                if (temp == 0) {
                                    p_node_neigh1 = p_node_neigh ;
                                    temp = 1 ;
                                }else{
                                    p_node_neigh2 = p_node_neigh ;
                                }  
                            }
                        }
                        auto& r_node1 = *p_node_neigh1 ; 
                        auto& r_node2 = *p_node_neigh2 ;
                        KRATOS_WATCH(r_node1)
                        KRATOS_WATCH(r_node2)
                        // Cerchiamo i vicini dei vicini
                        auto& r_node_neigh_vect_TOT = r_node1.GetValue(NEIGHBOUR_NODES);
                        auto& r_node_neigh_vect_2 = r_node2.GetValue(NEIGHBOUR_NODES);
                        KRATOS_WATCH(r_node_neigh_vect_TOT.size())
                        KRATOS_WATCH(r_node_neigh_vect_2.size()) 
                        // push_back
                        for(auto& it : r_node_neigh_vect_2.GetContainer()) {
                            r_node_neigh_vect_TOT.push_back(it);
                        }     
                        KRATOS_WATCH(r_node_neigh_vect_TOT.size())
                        // Unique
                        r_node_neigh_vect_TOT.Unique() ;    
                        KRATOS_WATCH(r_node_neigh_vect_TOT.size())
                        // shrink_to_fit
                        r_node_neigh_vect_TOT.shrink_to_fit() ;
                        KRATOS_WATCH(r_node_neigh_vect_TOT.size())
                        
                        // Inizio processo per calcolo gradiente
                        for (std::size_t i_neigh = 0; i_neigh < r_node_neigh_vect_TOT.size(); ++i_neigh) { 
                            auto p_node_neigh = r_node_neigh_vect_TOT(i_neigh).get();
                            double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                            KRATOS_WATCH(phi_distance)
                            if (phi_distance > 0 && p_node_neigh->IsNot(BOUNDARY) ) {
                                auto& r_node = *p_node_neigh ;
                                KRATOS_WATCH(r_node.Id())
                                double w_total = 0.0;
                                auto& r_elem_neigh_vect = r_node.GetValue(NEIGHBOUR_ELEMENTS);
                                // Set an auxilary map to calculate the current node nodal gradient contributions
                                std::map<std::size_t, Vector> neigh_dn_dx_map;
                                // Inizio calcolo coefficienti
                                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                                    if (p_elem_neigh->IsNot(BOUNDARY)) {
                                        if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                                            // Calculate the current element weight
                                            const auto& r_geom = p_elem_neigh->GetGeometry();
                                            const std::size_t n_nodes = r_geom.PointsNumber();
                                            // const double w = r_geom.DomainSize() / n_nodes;
                                            auto center_coord = r_geom.Center();
                                            const double w = 1/( (center_coord.X()-r_bd_node.X())*(center_coord.X()-r_bd_node.X()) + (center_coord.Y()-r_bd_node.Y())*(center_coord.Y()-r_bd_node.Y()) ) ;
                                            w_total += w;
                                            // Calculate the current element nodal gradient
                                            // Note that in here we assume that the gradient is constant within the element (simplex geometries)
                                            ShapeFunctionsGradientsType grad_N;
                                            r_geom.ShapeFunctionsIntegrationPointsGradients(grad_N, GeometryData::IntegrationMethod::GI_GAUSS_1);

                                            // Save the current element nodal weights
                                            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                                                std::size_t aux_id = r_geom[i_node].Id();
                                                const Vector r_node_dn_dx = row(grad_N[0], i_node);
                                                auto it_found = neigh_dn_dx_map.find(aux_id);
                                                if (it_found != neigh_dn_dx_map.end()) {
                                                    auto& r_val = it_found->second;
                                                    r_val += w * r_node_dn_dx;
                                                } else {
                                                    neigh_dn_dx_map.insert(std::make_pair(aux_id, w*r_node_dn_dx));
                                                }
                                            }
                                        }
                                    }
                                }
                                // Divide current coefficients by the total weight
                                for (auto& neigh_data : neigh_dn_dx_map) {
                                    auto& r_val = neigh_data.second;
                                    r_val /= w_total;
                                    KRATOS_WATCH(r_val)
                                }
                                // Update the surrogate nodes data
                                sur_bd_nodes_map.insert(std::make_pair(r_bd_node.Id(), neigh_dn_dx_map));
                            }
                        }
                        KRATOS_WATCH(r_bd_node.Id())
                        // if (r_bd_node.Id() == 21232) {
                        //     exit(0) ; 
                        // }
                        
                    //  STARDARD PROBLEMATIC Node --> ha 3 vicini, but just one not(BOUNDARY)
                    } else {          
                        for (std::size_t i_neigh = 0; i_neigh < r_node_neigh_vect.size(); ++i_neigh) {
                            auto p_node_neigh = r_node_neigh_vect(i_neigh).get();
                            auto& r_node = *p_node_neigh ;
                            KRATOS_WATCH(r_node.Id())
                            double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                            if (phi_distance > 0 && p_node_neigh->IsNot(BOUNDARY) ) {
                                // I get the "internal" node & // Find the node
                                auto& r_node = *p_node_neigh ;
                                KRATOS_WATCH(r_node.Id())
                                // Compute the gradient at this node, without using the elements that
                                // have at least one surrogate node, i.e. without using the INTERFACE elements
                                //_______________________________________________________________________________
                                // Calculate the nodal gradient coefficients in the MARKER elements neighbouring current node
                                double w_total = 0.0;
                                auto& r_elem_neigh_vect = r_node.GetValue(NEIGHBOUR_ELEMENTS);
                                // Set an auxilary map to calculate the current node nodal gradient contributions
                                std::map<std::size_t, Vector> neigh_dn_dx_map;
                                // In the case of concave skin_model_part can happend that zero of the elements around 
                                // the "internal" node is available for compute the gradient --> CONCAVE_PROBLEMATIC node
                                // Add a check
                                int count_concave = 0 ;
                                for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                                    auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                                    if (p_elem_neigh->IsNot(BOUNDARY)) {
                                        if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                                            count_concave = count_concave + 1 ;
                                            // Calculate the current element weight
                                            const auto& r_geom = p_elem_neigh->GetGeometry();
                                            const std::size_t n_nodes = r_geom.PointsNumber();
                                            // const double w = r_geom.DomainSize() / n_nodes;
                                            auto center_coord = r_geom.Center();
                                            const double w = 1/( (center_coord.X()-r_bd_node.X())*(center_coord.X()-r_bd_node.X()) + (center_coord.Y()-r_bd_node.Y())*(center_coord.Y()-r_bd_node.Y()) ) ;
                                            w_total += w;
                                            // Calculate the current element nodal gradient
                                            // Note that in here we assume that the gradient is constant within the element (simplex geometries)
                                            ShapeFunctionsGradientsType grad_N;
                                            r_geom.ShapeFunctionsIntegrationPointsGradients(grad_N, GeometryData::IntegrationMethod::GI_GAUSS_1);

                                            // Save the current element nodal weights
                                            for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                                                std::size_t aux_id = r_geom[i_node].Id();
                                                const Vector r_node_dn_dx = row(grad_N[0], i_node);
                                                auto it_found = neigh_dn_dx_map.find(aux_id);
                                                if (it_found != neigh_dn_dx_map.end()) {
                                                    auto& r_val = it_found->second;
                                                    r_val += w * r_node_dn_dx;
                                                } else {
                                                    neigh_dn_dx_map.insert(std::make_pair(aux_id, w*r_node_dn_dx));
                                                }
                                            }
                                        }
                                    }
                                }
                                // just if the node is not a CONCAVE_PROBLEMATIC node
                                if (count_concave !=0 ) {
                                    // Divide current coefficients by the total weight
                                    for (auto& neigh_data : neigh_dn_dx_map) {
                                        auto& r_val = neigh_data.second;
                                        r_val /= w_total;
                                    
                                    }
                                    // Update the surrogate nodes data
                                    sur_bd_nodes_map.insert(std::make_pair(r_bd_node.Id(), neigh_dn_dx_map));
                                }
                                // Check if there are CONCAVE_PROBLEMATIC nodes :
                                if (count_concave == 0) {
                                    KRATOS_WATCH("CONCAVE_PROBLEMATIC node")
                                    KRATOS_WATCH(r_bd_node.Id())
                                    // KRATOS_WATCH(r_node.Id())
                                    // Get the neighbours node of the previous candidate r_node
                                    auto& r_node_neigh_vect = r_node.GetValue(NEIGHBOUR_NODES);
                                    // Count how many available nodes there are
                                    int count_neigh_nodes = 0 ; 
                                    // Count the neighbours nodes that have DISTANCE > 0 & that are not BOUNDARY
                                    for (std::size_t i_neigh = 0; i_neigh < r_node_neigh_vect.size(); ++i_neigh) {
                                        auto p_node_neigh = r_node_neigh_vect(i_neigh).get();
                                        double &phi_distance = p_node_neigh->GetSolutionStepValue(DISTANCE);
                                        if (phi_distance > 0 && p_node_neigh->IsNot(BOUNDARY)) {
                                            count_neigh_nodes = count_neigh_nodes + 1 ;
                                            // Get the available node
                                            auto& r_node_concave = *p_node_neigh ;
                                            // KRATOS_WATCH(r_node_concave.Id())
                                            // Iniziamo a calcolare il gradiente
                                            double w_total = 0.0;
                                            auto& r_elem_neigh_vect = r_node_concave.GetValue(NEIGHBOUR_ELEMENTS);
                                            std::map<std::size_t, Vector> neigh_dn_dx_map;
                                            // Count available elements
                                            int count_avail_elements = 0 ;
                                            for (std::size_t i_neigh = 0; i_neigh < r_elem_neigh_vect.size(); ++i_neigh) {
                                                auto p_elem_neigh = r_elem_neigh_vect(i_neigh).get();
                                                if (p_elem_neigh->IsNot(BOUNDARY)) {
                                                    if (p_elem_neigh != nullptr && p_elem_neigh->Is(MARKER)) {
                                                        count_avail_elements = count_avail_elements + 1 ;
                                                        const auto& r_geom = p_elem_neigh->GetGeometry();
                                                        const std::size_t n_nodes = r_geom.PointsNumber();
                                                        // const double w = r_geom.DomainSize() / n_nodes;
                                                        auto center_coord = r_geom.Center();
                                                        const double w = 1/( (center_coord.X()-r_bd_node.X())*(center_coord.X()-r_bd_node.X()) + (center_coord.Y()-r_bd_node.Y())*(center_coord.Y()-r_bd_node.Y()) ) ;
                                                        w_total += w;
                                                        ShapeFunctionsGradientsType grad_N;
                                                        r_geom.ShapeFunctionsIntegrationPointsGradients(grad_N, GeometryData::IntegrationMethod::GI_GAUSS_1);
                                                        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                                                            std::size_t aux_id = r_geom[i_node].Id();
                                                            const Vector r_node_dn_dx = row(grad_N[0], i_node);
                                                            auto it_found = neigh_dn_dx_map.find(aux_id);
                                                            if (it_found != neigh_dn_dx_map.end()) {
                                                                auto& r_val = it_found->second;
                                                                r_val += w * r_node_dn_dx;
                                                            } else {
                                                                neigh_dn_dx_map.insert(std::make_pair(aux_id, w*r_node_dn_dx));
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            KRATOS_WATCH(count_avail_elements)
                                            if (count_avail_elements < 1) {
                                                KRATOS_WATCH("Nessun elemento vicino al candidato è utilizzabile")
                                                exit(0) ;
                                            }
                                            // Solo per CONCAVE_PROBLEMATIC node
                                            // Divide current coefficients by the total weight
                                            for (auto& neigh_data : neigh_dn_dx_map) {
                                                auto& r_val = neigh_data.second;
                                                r_val /= w_total;
                                            }
                                            // Update the surrogate nodes data
                                            sur_bd_nodes_map.insert(std::make_pair(r_bd_node.Id(), neigh_dn_dx_map));
                                        }
                                    }
                                    KRATOS_WATCH(count_neigh_nodes)
                                    if (count_neigh_nodes == 1) {
                                        KRATOS_WATCH("Tutto bene il concave problematic node")
                                    } else {
                                        if (count_neigh_nodes > 1) {
                                            KRATOS_WATCH("Ci sono due vicini disponibili, quale scelgo?")
                                        }
                                        KRATOS_WATCH("Bel problema... The Concave problematic is super problematic")
                                        exit(0) ;
                                    }
                                }

                                

                            }  
                        }      
                    }
                }
            }
        }
        return sur_bd_nodes_map;
    }

} // namespace Kratos.
