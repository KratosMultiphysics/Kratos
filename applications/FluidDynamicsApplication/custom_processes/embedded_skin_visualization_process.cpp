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
#include <unordered_map>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"

// Application includes
#include "embedded_skin_visualization_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    const ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart)
    : Process(), mrModelPart(rModelPart), mrVisualizationModelPart(rVisualizationModelPart) {
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    const ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart), mrVisualizationModelPart(rVisualizationModelPart) {

    Parameters default_parameters( R"(
    {
        "model_part_name"               : "default_model_part_name"
        "visualization_model_part_name" : "visualization_default_model_part_name"
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);

}

void EmbeddedSkinVisualizationProcess::ExecuteInitialize() {

    KRATOS_TRY;

    // Required variables check
    const auto& r_node = *mrModelPart.NodesBegin();
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

    KRATOS_CATCH("");
}

void EmbeddedSkinVisualizationProcess::ExecuteAfterOutputStep() {

    // Initialize the ids. for the new entries
    unsigned int node_id = mrModelPart.NumberOfNodes() + 1;
    unsigned int elem_id = mrModelPart.NumberOfElements() + 1;
    unsigned int cond_id = mrModelPart.NumberOfConditions() + 1;

    // Create the origin model part unordered_map
    typedef boost::unordered_map<std::vector<unsigned int>, Geometry<Node<3>>::Pointer, KeyHasher, KeyComparor> HashMapType;
    HashMapType origin_model_part_map;

    //#pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem)
    {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        // Get element geometry
        Geometry<Node<3>> &r_geometry = it_elem->GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();

        // Check if the element is split
        double n_pos (0.0), n_neg(0.0);
        Vector nodal_distances(n_nodes);
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
            const double dist = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
            nodal_distances[i_node] = dist;
            if (dist > 0.0)
                n_pos++;
            else
                n_neg++;
        }
        bool is_split = (n_pos > 0 && n_neg > 0) ? true : false;


        // If is split, save it in the visualization model part
        if (is_split){

            // Get the no-split geometry face local ids
            Vector face_n_nodes;
            Matrix face_loc_ids;
            r_geometry.NodesInFaces(face_loc_ids);
            r_geometry.NumberNodesInFaces(face_n_nodes);
            std::size_t n_faces = r_geometry.FacesNumber();

            std::vector<bool> add_face_vector;
            add_face_vector.resize(n_faces);
            for (unsigned int i_face = 0; i_face < n_faces; ++i_face){
                
                // Get the face global ids from the local ones
                unsigned int aux_count = 0;
                std::vector<unsigned int> face_nodes_ids;
                const unsigned int n_nodes_face = face_n_nodes[i_face];
                for (unsigned int i = 0; i < n_nodes_face; ++i){
                    const unsigned int local_id = face_loc_ids(aux_count++,i_face);
                    face_nodes_ids[i] = r_geometry[local_id].Id();
                }

                // Sort the global ids array
                std::sort(face_nodes_ids.begin(),face_nodes_ids.end());

                // Search for the current face in the origin model part map
                hashmap::iterator r_neighbour_geom = faces_map.find(face_nodes_ids);
                if (r_neighbour_geom != faces_map.end()){
                    add_face_vector[i_face] = false;
                } else {
                    add_face_vector[i_face] = true;
                }

                // If the face is already added, I know its global ids
                
                // If not, some nodes might be already added check that node by node
                
            }

            


            // Set the split utility and compute the splitting pattern
            DivideGeometry::Pointer p_split_utility = this->GetGeometrySplitUtility(r_geometry, nodal_distances);
            p_split_utility->GenerateDivision();
            p_split_utility->GenerateIntersectionsSkin();

            // Save the geometries from the splitting pattern in the visualization model part
            const unsigned int n_pos_geom = (p_split_utility->mPositiveSubdivisions).size();
            const unsigned int n_neg_geom = (p_split_utility->mNegativeSubdivisions).size();

            for (unsigned int i_geom = 0; i_geom < n_pos_geom; ++i_geom){
                auto p_geometry = (p_split_utility->mPositiveSubdivisions)[i_geom];
                const unsigned int n_nodes = p_geometry->PointsNumber();

                // Get subgeometry edges
                auto edges_vect = p_geometry->Edges();
                std::size_t n_edges = edges_vect.size(); 

                // Iterate the edges
                for (unsigned int i_edge = 0; i_edge < n_edges; ++i_edge){
                    auto p_edge = edges_vect[i_edge];
                }

    
                // Create the new nodes in the visualization model part
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
                    // Check if the node already exists
                    bool node_exists = false;
                    // TODO: IMPLEMENT IN HERE THE std::unordered_map CHECK

                    // If the node does not exist, it is created with constructor which takes a source node is used
                    if (!node_exists){
                        //mrVisualizationModelPart.CreateNewNode(node_id, p_geometry->GetPoint(i_node), 0);
                        // TODO: ADD THE NODE TO THE std::unordered_map
                        auto &r_orig_node = p_geometry->GetPoint(i_node);
                        Node<3>::Pointer p_new_node = mrVisualizationModelPart.CreateNewNode(node_id, r_orig_node.X(), r_orig_node.Y(), r_orig_node.Z());
                        node_id++;

                    }

                }

                // Create a new element with the subgeometry in the visualization model part
            }



        }



    }
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/
DivideGeometry::Pointer EmbeddedSkinVisualizationProcess::GetGeometrySplitUtility(
    const Geometry<Node<3>> &rGeometry,
    const Vector& rNodalDistances){

    // Get the geometry type
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();

    // Return the split utility according to the geometry tyoe
    switch (geometry_type){
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return Kratos::make_shared<DivideTriangle2D3>(rGeometry, rNodalDistances);
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return Kratos::make_shared<DivideTetrahedra3D4>(rGeometry, rNodalDistances);
        default:
            KRATOS_ERROR << "Asking for non-implemented geometry splitting utility.";
    }   
}

};  // namespace Kratos.
