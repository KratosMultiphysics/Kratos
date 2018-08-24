//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

// Application includes
#include "embedded_skin_visualization_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    const std::vector<Variable< double> >& rVisualizationScalarVariables,
    const std::vector<Variable< array_1d<double, 3> > >& rVisualizationVectorVariables,
    const std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >& rVisualizationComponentVariables,
    const std::string& rShapeFunctions,
    const bool ReformModelPartAtEachTimeStep) :
    Process(),
    mrModelPart(rModelPart),
    mrVisualizationModelPart(rVisualizationModelPart),
    mVisualizationScalarVariables(rVisualizationScalarVariables),
    mVisualizationVectorVariables(rVisualizationVectorVariables),
    mVisualizationComponentVariables(rVisualizationComponentVariables),
    mShapeFunctions(rShapeFunctions),
    mReformModelPartAtEachTimeStep(ReformModelPartAtEachTimeStep){
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    Parameters& rParameters) :
    Process(),
    mrModelPart(rModelPart),
    mrVisualizationModelPart(rVisualizationModelPart) {

    Parameters default_parameters( R"(
    {
        "shape_functions"                     : "standard",
        "reform_model_part_at_each_time_step" : false,
        "visualization_variables"             : ["VELOCITY","PRESSURE","DISTANCE"]
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

    // Validate if given shape functions value is admissible
    std::set<std::string> available_shape_functions = {"standard","ausas"};
    std::stringstream error_msg;

    if (available_shape_functions.find(rParameters["shape_functions"].GetString()) == available_shape_functions.end()){
        error_msg << "currently prescribed shape_functions : " << rParameters["shape_functions"].GetString() << std::endl;
        error_msg << "\tAdmissible values are : standard, ausas" << std::endl;
        KRATOS_ERROR << error_msg.str();
    }

    mShapeFunctions = rParameters["shape_functions"].GetString();
    mReformModelPartAtEachTimeStep = rParameters["reform_model_part_at_each_time_step"].GetBool();

    const unsigned int n_variables = rParameters["visualization_variables"].size();
    for (unsigned int i_var = 0; i_var < n_variables; ++i_var){
        std::string variable_name = rParameters["visualization_variables"].GetArrayItem(i_var).GetString();

        if(KratosComponents<Variable<double>>::Has(variable_name)){
            mVisualizationScalarVariables.push_back(KratosComponents< Variable<double> >::Get(variable_name));
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
            mVisualizationVectorVariables.push_back(KratosComponents< Variable< array_1d<double, 3> > >::Get(variable_name));
        } else if (KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >::Has(variable_name)) {
            mVisualizationComponentVariables.push_back(KratosComponents< VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >::Get(variable_name));
        } else {
            KRATOS_ERROR << "Only double, component and vector variables are allowed in the visualization variables list." ;
        }
    }
}

void EmbeddedSkinVisualizationProcess::ExecuteInitialize() {
    KRATOS_TRY;

    // Check that model part is not empty
    KRATOS_ERROR_IF(mrModelPart.Nodes().size() == 0) << "There are no nodes in the origin model part.";
    KRATOS_ERROR_IF(mrModelPart.Elements().size() == 0) << "There are no elements in the origin model part.";

    // Required variables check
    const auto &r_orig_node = *mrModelPart.NodesBegin();

    // Check visualization scalar variables
    for (unsigned int i_var = 0; i_var < mVisualizationScalarVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationScalarVariables[i_var], r_orig_node);
    }

    // Check visualization vector variables
    for (unsigned int i_var = 0; i_var < mVisualizationVectorVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationVectorVariables[i_var], r_orig_node);
    }

    // Check visualization component variables
    for (unsigned int i_var = 0; i_var < mVisualizationComponentVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationComponentVariables[i_var], r_orig_node);
    }

    KRATOS_CATCH("");
}

void EmbeddedSkinVisualizationProcess::ExecuteBeforeSolutionLoop() {
    if (mSetVisualizationMesh){
        // Copy the original nodes to the visualization model part
        this->CopyOriginNodes();

        // Creates the visualization model part geometrical entities (elements and conditions)
        this->CreateVisualizationGeometries();

        // Avoid creating the visualization mesh again
        mSetVisualizationMesh = false;
    }
}

void EmbeddedSkinVisualizationProcess::ExecuteInitializeSolutionStep() {
    if (mReformModelPartAtEachTimeStep){
        this->ExecuteBeforeSolutionLoop();
    }
}

void EmbeddedSkinVisualizationProcess::ExecuteBeforeOutputStep() {
    // Copy the origin model part nodal values for the non-intersection nodes
    this->CopyOriginNodalValues();

    // Compute new intersection nodes interpolation
    this->ComputeNewNodesInterpolation();
}

void EmbeddedSkinVisualizationProcess::ExecuteFinalizeSolutionStep(){
    if (mReformModelPartAtEachTimeStep){
        // Clear the new nodes hash map
        mCutNodesMap.clear();

        // Clear the new elements pointers vector
        mNewElementsPointers.clear();

        // Mark all the nodes, elements and conditions to be removed
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(TO_ERASE, true, mrVisualizationModelPart.Nodes());
        VariableUtils().SetFlag<ModelPart::ElementsContainerType>(TO_ERASE, true, mrVisualizationModelPart.Elements());
        VariableUtils().SetFlag<ModelPart::ConditionsContainerType>(TO_ERASE, true, mrVisualizationModelPart.Conditions());

        // Remove all the nodes, elements and conditions
        mrVisualizationModelPart.RemoveNodes(TO_ERASE);
        mrVisualizationModelPart.RemoveElements(TO_ERASE);
        mrVisualizationModelPart.RemoveConditions(TO_ERASE);

        // Initialize the create visualization mesh flag again
        mSetVisualizationMesh = true;
    }
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

void EmbeddedSkinVisualizationProcess::ComputeNewNodesInterpolation(){
    // For all the new elements, compute the interpolation with the proper shape functions
    // Note that this can be done in parallel since the intersection nodes are duplicated
    const int n_new_elems = mNewElementsPointers.size();

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < n_new_elems; ++i_elem){
        // Get element geometry
        auto it_elem = mNewElementsPointers.begin() + i_elem;
        const auto &r_geometry = it_elem->GetGeometry();
        const unsigned int n_points = r_geometry.PointsNumber();

        // For the generated elements, set the new interpolation values
        // Bear in mind that the non-split geometries kept the parent geometry
        // values, since they have been built using the origin model part nodes.
        for (unsigned int i_node = 0; i_node < n_points; ++i_node){
            // Search for the current node in the cut nodes hash map
            auto p_node = r_geometry(i_node);
            const CutNodesMapType::iterator cut_node_info = mCutNodesMap.find(p_node);

            // If it is not find, keep the origin value
            // Otherwise, the nodal values are computed with the modified shape functions values
            // Note that there is no necessity of locking the nodes since are always repeated
            if (cut_node_info != mCutNodesMap.end()){

                // Get the cut node info from the map tuple iterator
                Node<3>::Pointer p_edge_node_i, p_edge_node_j;
                double weight_edge_node_i, weight_edge_node_j;
                std::tie(p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j) = std::get<1>(*cut_node_info);

                // Interpolate the scalar variables
                for (unsigned int i_var = 0; i_var < mVisualizationScalarVariables.size(); ++i_var){
                    const Variable<double> &r_scalar_var = mVisualizationScalarVariables[i_var];
                    const double &edge_node_i_value = p_edge_node_i->FastGetSolutionStepValue(r_scalar_var);
                    const double &edge_node_j_value = p_edge_node_j->FastGetSolutionStepValue(r_scalar_var);
                    p_node->FastGetSolutionStepValue(r_scalar_var) = weight_edge_node_i * edge_node_i_value + weight_edge_node_j * edge_node_j_value;
                }

                // Interpolate the vector variables
                for (unsigned int i_var = 0; i_var < mVisualizationVectorVariables.size(); ++i_var){
                    const Variable<array_1d<double, 3> > &r_vector_var = mVisualizationVectorVariables[i_var];
                    const array_1d<double, 3> &edge_node_i_value = p_edge_node_i->FastGetSolutionStepValue(r_vector_var);
                    const array_1d<double, 3> &edge_node_j_value = p_edge_node_j->FastGetSolutionStepValue(r_vector_var);
                    p_node->FastGetSolutionStepValue(r_vector_var) = weight_edge_node_i * edge_node_i_value + weight_edge_node_j * edge_node_j_value;
                }

                // Interpolate the component variables
                for (unsigned int i_var = 0; i_var < mVisualizationComponentVariables.size(); ++i_var){
                    const VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > &r_comp_var = mVisualizationComponentVariables[i_var];
                    const double &edge_node_i_value = p_edge_node_i->FastGetSolutionStepValue(r_comp_var);
                    const double &edge_node_j_value = p_edge_node_j->FastGetSolutionStepValue(r_comp_var);
                    p_node->FastGetSolutionStepValue(r_comp_var) = weight_edge_node_i * edge_node_i_value + weight_edge_node_j * edge_node_j_value;
                }
            }
        }
    }
}

void EmbeddedSkinVisualizationProcess::CopyOriginNodes(){
    // Creates a copy of all the origin model part nodes to the visualization model part
    // Note that the original nodes will be reused when creating the splitting geometries
    const int n_nodes = mrModelPart.NumberOfNodes();
    ModelPart::NodeIterator orig_nodes_begin = mrModelPart.NodesBegin();
    for (int i_node = 0; i_node < n_nodes; ++i_node){
        auto it_node = orig_nodes_begin + i_node;
        auto p_node = mrVisualizationModelPart.CreateNewNode(it_node->Id(), *it_node);
    }
}

void EmbeddedSkinVisualizationProcess::CopyOriginNodalValues(){
    const unsigned int n_old_nodes = mrModelPart.NumberOfNodes();

    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(n_old_nodes); ++i_node){
        const auto it_old_node = mrModelPart.NodesBegin() + i_node;
        auto it_new_node = mrVisualizationModelPart.NodesBegin() + i_node;

        for (unsigned int i_var = 0; i_var < mVisualizationScalarVariables.size(); ++i_var){
            it_new_node->FastGetSolutionStepValue(mVisualizationScalarVariables[i_var]) =
                it_old_node->FastGetSolutionStepValue(mVisualizationScalarVariables[i_var]);
        }

        for (unsigned int i_var = 0; i_var < mVisualizationVectorVariables.size(); ++i_var){
            it_new_node->FastGetSolutionStepValue(mVisualizationVectorVariables[i_var]) =
                it_old_node->FastGetSolutionStepValue(mVisualizationVectorVariables[i_var]);
        }

        for (unsigned int i_var = 0; i_var < mVisualizationComponentVariables.size(); ++i_var){
            it_new_node->FastGetSolutionStepValue(mVisualizationComponentVariables[i_var]) =
                it_old_node->FastGetSolutionStepValue(mVisualizationComponentVariables[i_var]);
        }
    }
}

void EmbeddedSkinVisualizationProcess::CreateVisualizationGeometries(){

    int n_nodes = mrModelPart.NumberOfNodes();
    int n_elems = mrModelPart.NumberOfElements();
    int n_conds = mrModelPart.NumberOfConditions();

    // Initialize the ids. for the new entries
    // For the temporal ids. there is no necessity of synchronizing between processors
    unsigned int temp_node_id = (n_nodes > 0) ? ((mrModelPart.NodesEnd()-1)->Id()) + 1 : 1;
    unsigned int temp_elem_id = (n_elems > 0) ? ((mrModelPart.ElementsEnd()-1)->Id()) + 1 : 1;
    unsigned int temp_cond_id = (n_conds > 0) ? ((mrModelPart.ConditionsEnd()-1)->Id()) + 1 : 1;

    // Auxiliar vectors to store pointers to the new entities
    // These vectors will be use when renumbering the entities ids.
    ModelPart::NodesContainerType new_nodes_vect;
    ModelPart::ConditionsContainerType new_conds_vect;

    // Set the new entities properties
    Properties::Pointer p_pos_prop, p_neg_prop;
    std::tie(p_pos_prop, p_neg_prop) = this->SetVisualizationProperties();

    // Add the elements to the visualization model part
    for (int i_elem = 0; i_elem < n_elems; ++i_elem){
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i_elem;

        // Get element geometry
        const Geometry<Node<3>>::Pointer p_geometry = it_elem->pGetGeometry();
        const unsigned int n_nodes = p_geometry->PointsNumber();
        const Vector nodal_distances = this->SetDistancesVector(it_elem);

        // Check if the element is split
        const bool is_split = this->ElementIsSplit(p_geometry, nodal_distances);

        // If the element is split, create the new entities
        if (is_split){
            // Set the split utility and compute the splitting pattern
            ModifiedShapeFunctions::Pointer p_modified_shape_functions = this->SetModifiedShapeFunctionsUtility(p_geometry, nodal_distances);
            DivideGeometry::Pointer p_split_utility = p_modified_shape_functions->pGetSplittingUtil();

            // Create the auxiliar map that will be used to generate the skin
            std::unordered_map<std::pair<unsigned int,bool>, unsigned int, Hash, KeyEqual> new_nodes_map;

            // Get the split geometries from the splitting pattern
            const unsigned int n_pos_split_geom = (p_split_utility->mPositiveSubdivisions).size();
            const unsigned int n_neg_split_geom = (p_split_utility->mNegativeSubdivisions).size();
            std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_geometries;
            split_geometries.reserve(n_pos_split_geom + n_neg_split_geom);
            split_geometries.insert(split_geometries.end(), (p_split_utility->mPositiveSubdivisions).begin(), (p_split_utility->mPositiveSubdivisions).end());
            split_geometries.insert(split_geometries.end(), (p_split_utility->mNegativeSubdivisions).begin(), (p_split_utility->mNegativeSubdivisions).end());

            // Create the split geometries in the visualization model part
            for (unsigned int i_geom = 0; i_geom < split_geometries.size(); ++i_geom){
                const bool pos_side = i_geom < n_pos_split_geom ? true : false;
                const DivideGeometry::IndexedPointGeometryPointerType p_sub_geom = split_geometries[i_geom];
                const unsigned int sub_geom_n_nodes = p_sub_geom->PointsNumber();

                // Fill the new element nodes array
                Element::NodesArrayType sub_geom_nodes_array;
                for (unsigned int i_sub_geom_node = 0; i_sub_geom_node < sub_geom_n_nodes; ++i_sub_geom_node){

                    DivideGeometry::IndexedPointType &sub_geom_node = p_sub_geom->operator[](i_sub_geom_node);
                    const unsigned int local_id = sub_geom_node.Id();

                    // Non-intersection node. Get the copied original geometry nodes
                    // Note that we do not reuse the original model part node since this forces
                    // the new intersection nodes to have the same buffer size, variables list, etc.
                    if (local_id < sub_geom_n_nodes){
                        const unsigned int aux_gl_id = (p_geometry->operator()(local_id))->Id();
                        sub_geom_nodes_array.push_back(mrVisualizationModelPart.pGetNode(aux_gl_id));
                        // Intersection node. The node is located in an intersected edge.
                        // Thus, take it from the geometry created by the splitting util.
                    } else {
                        const unsigned int intersected_edge_id = local_id - n_nodes;

                        // Get the intersected edge node_i and node_j
                        const unsigned int node_i = (p_split_utility->GetEdgeIdsI())[intersected_edge_id];
                        const unsigned int node_j = (p_split_utility->GetEdgeIdsJ())[intersected_edge_id];

                        // Create a new node based in the indexed point auxiliar node
                        const array_1d<double, 3> point_coords = sub_geom_node.Coordinates();
                        Node<3>::Pointer p_new_node = mrVisualizationModelPart.CreateNewNode(temp_node_id, point_coords[0], point_coords[1], point_coords[2]);
                        sub_geom_nodes_array.push_back(p_new_node);
                        new_nodes_vect.push_back(p_new_node);

                        // Add the new node info to the hash map
                        const Node<3>::Pointer p_node_i = p_geometry->operator()(node_i);
                        const Node<3>::Pointer p_node_j = p_geometry->operator()(node_j);

                        Matrix edge_N_values;
                        if (i_geom < n_pos_split_geom){
                            p_modified_shape_functions->ComputeShapeFunctionsOnPositiveEdgeIntersections(edge_N_values);
                        } else {
                            p_modified_shape_functions->ComputeShapeFunctionsOnNegativeEdgeIntersections(edge_N_values);
                        }

                        const double node_i_weight = edge_N_values(intersected_edge_id, node_i);
                        const double node_j_weight = edge_N_values(intersected_edge_id, node_j);

                        auto new_node_info = std::make_tuple(p_node_i, p_node_j, node_i_weight, node_j_weight);
                        mCutNodesMap.insert(CutNodesMapType::value_type(p_new_node, new_node_info));

                        // Link the new node global id. to its local index in the splitting util. Note that the side
                        // (positive or negative) is saved as well. This will be used when setting up the skin conditions.
                        std::pair<unsigned int, bool> aux_info(local_id,pos_side);
                        std::pair<std::pair<unsigned int,bool>, unsigned int> new_pair(aux_info, temp_node_id);
                        new_nodes_map.insert(new_pair);

                        // Update the new nodes id. counter
                        temp_node_id++;
                    }
                }

                // Set the new element properties
                Properties::Pointer p_elem_prop = pos_side ? p_pos_prop : p_neg_prop;

                // Create the new element
                Element::Pointer p_new_elem = it_elem->Create(temp_elem_id, sub_geom_nodes_array, p_elem_prop);
                mNewElementsPointers.push_back(p_new_elem);

                // Update the new elements id. counter
                temp_elem_id++;
            }

            // Get the interface geometries from the splitting pattern
            const unsigned int n_pos_interface_geom = (p_split_utility->mPositiveInterfaces).size();
            const unsigned int n_neg_interface_geom = (p_split_utility->mNegativeInterfaces).size();

            std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_interface_geometries;
            split_interface_geometries.reserve(n_pos_interface_geom + n_neg_interface_geom);
            split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mPositiveInterfaces).begin(), (p_split_utility->mPositiveInterfaces).end());
            split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mNegativeInterfaces).begin(), (p_split_utility->mNegativeInterfaces).end());

            // Create the split interface geometries in the visualization model part
            for (unsigned int i_int_geom = 0; i_int_geom < split_interface_geometries.size(); ++i_int_geom){
                const bool int_pos_side = (i_int_geom < n_pos_interface_geom) ? true : false;
                DivideGeometry::IndexedPointGeometryPointerType p_int_sub_geom = split_interface_geometries[i_int_geom];
                GeometryData::KratosGeometryType p_int_sub_geom_type = p_int_sub_geom->GetGeometryType();
                const unsigned int sub_int_geom_n_nodes = p_int_sub_geom->PointsNumber();

                // Fill the new condition nodes array
                Condition::NodesArrayType sub_int_geom_nodes_array;
                for (unsigned int i_node = 0; i_node < sub_int_geom_n_nodes; ++i_node){

                    DivideGeometry::IndexedPointType &sub_int_geom_node = p_int_sub_geom->operator[](i_node);
                    const unsigned int local_id = sub_int_geom_node.Id();

                    // Get the global id from the intersection nodes map
                    unsigned int global_id;
                    std::pair<unsigned int,bool> aux_int_info(local_id,int_pos_side);
                    auto got = new_nodes_map.find(aux_int_info);
                    if (got != new_nodes_map.end()){
                        global_id = got->second;
                    } else {
                        const std::string side = int_pos_side ? "positive" : "negative";
                        KRATOS_ERROR << "Local id " << std::get<0>(aux_int_info) << " in " << side << " side not found in new nodes map for element " << it_elem->Id();
                    }

                    sub_int_geom_nodes_array.push_back(mrVisualizationModelPart.pGetNode(global_id));
                }

                // Set the new condition geometry
                Geometry< Node<3> >::Pointer p_new_geom = SetNewConditionGeometry(
                    p_int_sub_geom_type,
                    sub_int_geom_nodes_array);

                // Set the new condition properties
                Properties::Pointer p_cond_prop = (i_int_geom < n_pos_interface_geom)? p_pos_prop : p_neg_prop;

                // Create the new condition
                Condition::Pointer p_new_cond = Kratos::make_shared<Condition>(temp_cond_id, p_new_geom, p_cond_prop);
                new_conds_vect.push_back(p_new_cond);
                mrVisualizationModelPart.AddCondition(p_new_cond);

                // Update the new elements id. counter
                temp_cond_id++;

            }
        // Otherwise add an element with the original geometry
        } else {
            // Copy the current element as it was in the origin model part according to its distance sign
            const bool is_positive = this->ElementIsPositive(p_geometry, nodal_distances);
            if (is_positive){
                mrVisualizationModelPart.AddElement(it_elem->Create(it_elem->Id(), p_geometry, p_pos_prop));
            } else {
                mrVisualizationModelPart.AddElement(it_elem->Create(it_elem->Id(), p_geometry, p_neg_prop));
            }
        }
    }

    // Once all the entities have been created, renumber the ids.
    // Created entities local number partial reduction
    int n_nodes_local = new_nodes_vect.size();
    int n_conds_local = new_conds_vect.size();
    int n_elems_local = mNewElementsPointers.size();
    int n_nodes_local_scansum, n_elems_local_scansum, n_conds_local_scansum;
    mrModelPart.GetCommunicator().ScanSum(n_nodes_local, n_nodes_local_scansum);
    mrModelPart.GetCommunicator().ScanSum(n_elems_local, n_elems_local_scansum);
    mrModelPart.GetCommunicator().ScanSum(n_conds_local, n_conds_local_scansum);

    // Origin model part number of entities
    int n_nodes_orig = mrModelPart.NumberOfNodes();
    int n_elems_orig = mrModelPart.NumberOfElements();
    int n_conds_orig = mrModelPart.NumberOfConditions();
    mrModelPart.GetCommunicator().SumAll(n_nodes_orig);
    mrModelPart.GetCommunicator().SumAll(n_elems_orig);
    mrModelPart.GetCommunicator().SumAll(n_conds_orig);

    // Initialize the new ids. values
    std::size_t new_node_id(n_nodes_orig + n_nodes_local_scansum - n_nodes_local + 1);
    std::size_t new_elem_id(n_elems_orig + n_elems_local_scansum - n_elems_local + 1);
    std::size_t new_cond_id(n_conds_orig + n_conds_local_scansum - n_conds_local + 1);

    // Set the new entities ids. values
    auto new_nodes_begin = new_nodes_vect.begin();
    auto new_conds_begin = new_conds_vect.begin();
    auto new_elems_begin = mNewElementsPointers.begin();

    for (int i_node = 0; i_node < static_cast<int>(new_nodes_vect.size()); ++i_node){
        auto it_node = new_nodes_begin + i_node;
        const unsigned int new_id = new_node_id + i_node;
        it_node->SetId(new_id);
    }

    for (int i_cond = 0; i_cond < static_cast<int>(new_conds_vect.size()); ++i_cond){
        auto it_cond = new_conds_begin + i_cond;
        const unsigned int new_id = new_cond_id + i_cond;
        it_cond->SetId(new_id);
    }

    for (int i_elem = 0; i_elem < static_cast<int>(mNewElementsPointers.size()); ++i_elem){
        auto it_elem = new_elems_begin + i_elem;
        const unsigned int new_id = new_elem_id + i_elem;
        it_elem->SetId(new_id);
    }

    mrVisualizationModelPart.AddElements(mNewElementsPointers.begin(), mNewElementsPointers.end());
    mrVisualizationModelPart.AddConditions(new_conds_vect.begin(), new_conds_vect.end());

    // Wait for all nodes to renumber its nodes
    mrModelPart.GetCommunicator().Barrier();
}

bool EmbeddedSkinVisualizationProcess::ElementIsPositive(
    Geometry<Node<3>>::Pointer pGeometry,
    const Vector &rNodalDistances){

    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        if (rNodalDistances[i_node] > 0.0)
            n_pos++;
    }

    const bool is_positive = (n_pos == pts_number) ? true : false;

    return is_positive;
}

bool EmbeddedSkinVisualizationProcess::ElementIsSplit(
    Geometry<Node<3>>::Pointer pGeometry,
    const Vector &rNodalDistances){

    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0), n_neg(0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        if (rNodalDistances[i_node] > 0.0)
            n_pos++;
        else
            n_neg++;
    }

    const bool is_split = (n_pos > 0 && n_neg > 0) ? true : false;

    return is_split;
}

const Vector EmbeddedSkinVisualizationProcess::SetDistancesVector(ModelPart::ElementIterator ItElem){

    auto &r_geom = ItElem->GetGeometry();
    Vector nodal_distances(r_geom.PointsNumber());

    if (mShapeFunctions == "standard"){
        // Continuous nodal distance function case
        for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            nodal_distances[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
        }
    } else if (mShapeFunctions == "ausas") {
        // Discontinuous elemental distance function case
        nodal_distances = ItElem->GetValue(ELEMENTAL_DISTANCES);
    } else {
        KRATOS_ERROR << "Asking for a non-implemented modified shape functions type.";
    }

    return nodal_distances;
}

ModifiedShapeFunctions::Pointer EmbeddedSkinVisualizationProcess::SetModifiedShapeFunctionsUtility(
    const Geometry<Node<3>>::Pointer pGeometry,
    const Vector& rNodalDistances){

    // Get the geometry type
    const GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();

    // Return the modified shape functions utility
    if (mShapeFunctions == "standard"){
        switch (geometry_type){
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);
            default:
                KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
        }
    } else if (mShapeFunctions == "ausas"){
        switch (geometry_type){
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                return Kratos::make_shared<Triangle2D3AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                return Kratos::make_shared<Tetrahedra3D4AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
            default:
                KRATOS_ERROR << "Asking for a non-implemented Ausas modified shape functions geometry.";
        }
    } else {
        KRATOS_ERROR << "Asking for a non-implemented modified shape functions type.";
    }
}

Geometry< Node<3> >::Pointer EmbeddedSkinVisualizationProcess::SetNewConditionGeometry(
    const GeometryData::KratosGeometryType &rOriginGeometryType,
    const Condition::NodesArrayType &rNewNodesArray){

    switch(rOriginGeometryType){
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            return Kratos::make_shared<Line2D2< Node<3> > >(rNewNodesArray);
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return Kratos::make_shared<Triangle3D3< Node<3> > >(rNewNodesArray);
        default:
            KRATOS_ERROR << "Implement the visualization for the intersection geometry type " << rOriginGeometryType;
    }
}

std::tuple< Properties::Pointer , Properties::Pointer > EmbeddedSkinVisualizationProcess::SetVisualizationProperties(){
    // Set the properties for the new elements depending if the
    // element is in the positive or negative side of the cut.
    // In this way, two layers will appear in GiD.
    unsigned int max_prop_id = 0;
    for (auto it_prop = mrModelPart.GetRootModelPart().PropertiesBegin(); it_prop < mrModelPart.GetRootModelPart().PropertiesEnd(); ++it_prop){
        if (max_prop_id < it_prop->Id()){
            max_prop_id = it_prop->Id();
        }
    }
    Properties::Pointer p_pos_prop = Kratos::make_shared<Properties>(max_prop_id + 1);
    Properties::Pointer p_neg_prop = Kratos::make_shared<Properties>(max_prop_id + 2);
    mrVisualizationModelPart.AddProperties(p_pos_prop);
    mrVisualizationModelPart.AddProperties(p_neg_prop);

    return std::make_tuple(p_pos_prop , p_neg_prop);
}

};  // namespace Kratos.
