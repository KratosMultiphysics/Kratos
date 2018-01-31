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
    const std::vector<Variable< double> > VisualizationScalarVariables,
    const std::vector<Variable< array_1d<double, 3> > > VisualizationVectorVariables,
    const std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > > VisualizationComponentVariables,
    const std::string ShapeFunctions) : 
    Process(), 
    mrModelPart(rModelPart), 
    mrVisualizationModelPart(rVisualizationModelPart) , 
    mVisualizationScalarVariables(VisualizationScalarVariables),
    mVisualizationVectorVariables(VisualizationVectorVariables),
    mVisualizationComponentVariables(VisualizationComponentVariables),
    mShapeFunctions(ShapeFunctions){ 
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
        "shape_functions"         : "standard",
        "visualization_variables" : ["VELOCITY","PRESSURE"]
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mShapeFunctions = rParameters["shape_functions"].GetString();

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

    // Required variables check
    const auto& r_node = *mrModelPart.NodesBegin();
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

    // Check visualization scalar variables
    for (unsigned int i_var = 0; i_var < mVisualizationScalarVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationScalarVariables[i_var], r_node);
    }

    // Check visualization vector variables
    for (unsigned int i_var = 0; i_var < mVisualizationVectorVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationVectorVariables[i_var], r_node);
    }

    // Check visualization component variables
    for (unsigned int i_var = 0; i_var < mVisualizationComponentVariables.size(); ++i_var){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVisualizationComponentVariables[i_var], r_node);
    }

    KRATOS_CATCH("");
}

void EmbeddedSkinVisualizationProcess::ExecuteInitializeSolutionStep() {

    if (!mMeshCreationWasPerformed){

        const unsigned int n_nodes = mrModelPart.NumberOfNodes();
        const unsigned int n_elems = mrModelPart.NumberOfElements();
        const unsigned int n_conds = mrModelPart.NumberOfConditions();

        // Initialize the ids. for the new entries
        unsigned int new_node_id = n_nodes + 1;
        unsigned int new_elem_id = n_elems + 1;
        unsigned int new_cond_id = n_conds + 1;

        // Set the properties for the new elements depending if the 
        // element is in the positive or negative side of the cut.
        // In this way, two layers will appear in GiD.
        const unsigned int last_prop = ((*(mrModelPart.pProperties())).end()-1)->Id();
        Properties::Pointer p_pos_prop = Kratos::make_shared<Properties>(last_prop + 1);
        Properties::Pointer p_neg_prop = Kratos::make_shared<Properties>(last_prop + 2);
        mrVisualizationModelPart.AddProperties(p_pos_prop);
        mrVisualizationModelPart.AddProperties(p_neg_prop);

        // Add all the origin model part nodes to the visualization model part
        // Note that the original nodes will be reused when creating the splitting geometries
        mrVisualizationModelPart.AddNodes(mrModelPart.NodesBegin(), mrModelPart.NodesEnd());

        // Add the elements to the visualization model part
        auto elems_begin = mrModelPart.ElementsBegin();
        for (unsigned int i_elem = 0; i_elem < n_elems; ++i_elem){
            auto it_elem = elems_begin + i_elem;

            // Get element geometry
            Geometry<Node<3>>::Pointer p_geometry = it_elem->pGetGeometry();
            const unsigned int n_nodes = p_geometry->PointsNumber();
            Vector nodal_distances(n_nodes);

            // Check if the element is split
            const bool is_split = this->ElementIsSplit(p_geometry, nodal_distances);

            // If the element is split, call the splitting utility 
            if (is_split){
                // Set the split utility and compute the splitting pattern
                ModifiedShapeFunctions::Pointer p_modified_shape_functions = this->SetModifiedShapeFunctionsUtility(p_geometry, nodal_distances);
                DivideGeometry::Pointer p_split_utility = p_modified_shape_functions->pGetSplittingUtil();

                // Save the geometries from the splitting pattern in the visualization model part
                const unsigned int n_pos_split_geom = (p_split_utility->mPositiveSubdivisions).size();
                std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_geometries;
                split_geometries.insert(split_geometries.end(), (p_split_utility->mPositiveSubdivisions).begin(), (p_split_utility->mPositiveSubdivisions).end());
                split_geometries.insert(split_geometries.end(), (p_split_utility->mNegativeSubdivisions).begin(), (p_split_utility->mNegativeSubdivisions).end());

                for (unsigned int i_geom = 0; i_geom < split_geometries.size(); ++i_geom){
                    auto p_sub_geom = split_geometries[i_geom];
                    const unsigned int sub_geom_n_nodes = p_sub_geom->PointsNumber();

                    // Fill the new element nodes array
                    Element::NodesArrayType sub_geom_nodes_array;
                    for (unsigned int i_sub_geom_node = 0; i_sub_geom_node < sub_geom_n_nodes; ++i_sub_geom_node){

                        auto sub_geom_node = p_sub_geom->operator[](i_sub_geom_node);
                        const unsigned int local_id = sub_geom_node.Id();

                        // Non-intersection node. Take the node from the parent geometry
                        if (local_id < sub_geom_n_nodes){
                            sub_geom_nodes_array.push_back(p_geometry->operator()(local_id));

                        // Intersection node. The node is located in an intersected edge.
                        // Thus, take it from the geometry created by the splitting util.
                        } else {
                            const unsigned int intersected_edge_id = local_id - n_nodes;

                            // Get the intersected edge node_i and node_j
                            const unsigned int node_i = (p_split_utility->GetEdgeIdsI())[intersected_edge_id];
                            const unsigned int node_j = (p_split_utility->GetEdgeIdsJ())[intersected_edge_id];

                            // Create a new node based in the indexed point auxiliar node
                            const array_1d<double, 3> point_coords = sub_geom_node.Coordinates();
                            Node<3>::Pointer p_new_node = mrVisualizationModelPart.CreateNewNode(new_node_id, point_coords[0], point_coords[1], point_coords[2]);
                            sub_geom_nodes_array.push_back(p_new_node);

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

                            // Update the new nodes id. counter
                            new_node_id++;
                        }
                    }

                    // Set the new element properties
                    Properties::Pointer p_elem_prop = (i_geom < n_pos_split_geom)? p_pos_prop : p_neg_prop;

                    // Create the new element
                    auto p_new_elem = it_elem->Create(new_elem_id, sub_geom_nodes_array, p_elem_prop);
                    mrVisualizationModelPart.AddElement(p_new_elem);
                    mNewElementsIds.push_back(new_elem_id);

                    // Update the new elements id. counter
                    new_elem_id++;
                }
                
            // Otherwise add an element with the original geometry
            } else {
                // Add the current element as it was in the origin model part
                mrVisualizationModelPart.AddElement(*it_elem.base());

                // Once the element has been added to the visualization model part
                // modify its properties according to its distance sign.
                const bool is_positive = this->ElementIsPositive(p_geometry);

                if (is_positive){
                    (mrVisualizationModelPart.ElementsEnd()-1)->SetProperties(p_pos_prop);
                } else {
                    (mrVisualizationModelPart.ElementsEnd()-1)->SetProperties(p_neg_prop);
                }
            }
        }

        mMeshCreationWasPerformed = true;
    }
}

void EmbeddedSkinVisualizationProcess::ExecuteBeforeOutputStep() {


    // For all the new elements, compute the interpolation with the proper shape functions
    const int n_new_elems = mNewElementsIds.size();

    // Sort the elements pointer vector set. Otherwise,the access is not 
    // threadsafe since it tries to shrink to fit the capacity when first accessing.
    mrVisualizationModelPart.Elements().Sort();

    #pragma omp parallel for firstprivate(n_new_elems)
    for (int i_elem = 0; i_elem < n_new_elems; ++i_elem){
        // Get element geometry
        const int elem_id = mNewElementsIds[i_elem];
        Element::Pointer p_new_elem = mrVisualizationModelPart.Elements()(elem_id);
        Geometry<Node<3>> &r_geometry = p_new_elem->GetGeometry();
        const unsigned int n_points = r_geometry.PointsNumber();

        // For the generated elements, set the new interpolation values
        // Bear in mind that the non-split geometries kept the parent geometry
        // values, since they have been built using the origin model part nodes.
        for (unsigned int i_node = 0; i_node < n_points; ++i_node){
            auto p_node = r_geometry(i_node);

            // Search for the current node in the cut nodes hash map
            CutNodesMapType::iterator cut_node_info = mCutNodesMap.find(p_node);

            // If it is not find, keep the origin value
            // Otherwise, the nodal values are computed with the modified shape functions values
            // Note that there is no necessity of locking the nodes since are always repeated
            if (cut_node_info != mCutNodesMap.end()){
                // Get the cut node info from the map tuple iterator
                Node<3>::Pointer p_edge_node_i, p_edge_node_j;
                double weight_edge_node_i, weight_edge_node_j;
                std::tie(p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j) = std::get<1>(*cut_node_info);

                // Interpolate the nodal values
                double &p_cut_node = p_node->FastGetSolutionStepValue(PRESSURE);
                const double &p_node_i = p_edge_node_i->FastGetSolutionStepValue(PRESSURE);
                const double &p_node_j = p_edge_node_j->FastGetSolutionStepValue(PRESSURE);
                p_cut_node = weight_edge_node_i*p_node_i + weight_edge_node_j*p_node_j;

                array_1d<double, 3> &v_cut_node = p_node->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3>& v_node_i = p_edge_node_i->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3>& v_node_j = p_edge_node_j->FastGetSolutionStepValue(VELOCITY);
                v_cut_node = weight_edge_node_i*v_node_i + weight_edge_node_j*v_node_j;
            }
        }
    }
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

const bool EmbeddedSkinVisualizationProcess::ElementIsPositive(
    Geometry<Node<3>>::Pointer pGeometry){

    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0), n_neg(0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        const double dist = (pGeometry->operator[](i_node)).FastGetSolutionStepValue(DISTANCE);
        if (dist > 0.0)
            n_pos++;
        else
            n_neg++;
    }

    const bool is_positive = (n_pos == pts_number) ? true : false;

    return is_positive;

}

const bool EmbeddedSkinVisualizationProcess::ElementIsSplit(
    Geometry<Node<3>>::Pointer pGeometry,
    Vector &rNodalDistances){

    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0), n_neg(0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        const double dist = (pGeometry->operator[](i_node)).FastGetSolutionStepValue(DISTANCE);
        rNodalDistances[i_node] = dist;
        if (dist > 0.0)
            n_pos++;
        else
            n_neg++;
    }

    const bool is_split = (n_pos > 0 && n_neg > 0) ? true : false;

    return is_split;

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

};  // namespace Kratos.
