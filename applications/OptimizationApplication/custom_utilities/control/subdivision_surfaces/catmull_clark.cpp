//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Bastian Devresse
//


// System includes
#include <numeric>

// Project includes
#include "expression/literal_flat_expression.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "containers/model.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "spatial_containers/bucket.h"
#include "spatial_containers/kd_tree.h"
#include "spaces/ublas_space.h"
#include "opensubdiv_utilities.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/intersection_utilities.h"
#include "input_output/vtk_output.h"

/// BSpline surfaces
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
#include "../applications/ShapeOptimizationApplication/custom_utilities/geometry_utilities.h"   // TODO include/transfer geometry_utilities from SOA to OptApp

// External includes
#include <fcpw/fcpw.h>

// Include base h
#include "catmull_clark.h"


namespace Kratos
{
// namespace SDSUtils
// {

// using IndexType = std::size_t;
// using NodeType = Node::NodeType;
using GeometryType = Geometry<Node>;
using ModelPart = ModelPart;
// using MeshType = ModelPart::MeshType;
using CoordinatesArrayType = Point::CoordinatesArrayType;
typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

// typedefs for nearest point search
// typedef Node NodeType;
typedef NodeType::Pointer NodeTypePointer;
typedef std::vector<NodeType::Pointer> NodeVector;
typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
typedef std::vector<double>::iterator DoubleVectorIterator;
typedef Kratos::Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
typedef Tree< KDTreePartition<BucketType> > KDTree;

typedef OpenSubdiv::v3_6_0::Far::ConstIndexArray ConstIndexArray;
typedef OpenSubdiv::Far::TopologyDescriptor Descriptor;

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

using SparseMatrixType = SparseSpaceType::MatrixType;

CatmullClarkSDS::CatmullClarkSDS(ModelPart& rControlPolygon, const ModelPart& rControlledMesh) : mrControlPolygon{rControlPolygon}, mrControlledMesh(rControlledMesh)
{}

CatmullClarkSDS::~CatmullClarkSDS() {}

IndexType GetOppositeIndex(ConstIndexArray ValueArray, const IndexType Value)
{
    IndexType index = ValueArray.FindIndex(Value);
    SizeType max_index = ValueArray.size() - 1;
    if (index + 2 > max_index) {
        return ValueArray[index - 2];
    }
    else return ValueArray[index + 2];
}

std::vector<IndexType> ReorderIndicesToStartAtValue(const IndexType Value, ConstIndexArray ValueArray)
{
    SizeType starting_index = ValueArray.FindIndex(Value);
    SizeType max_index = ValueArray.size() - 1;
    std::vector<IndexType> ReturnArray;

    for (SizeType i = 0; i <= max_index - starting_index; ++i) {
        ReturnArray.push_back( ValueArray[starting_index + i] );
    }
    for (SizeType i = 0; i < starting_index; ++i) {
        ReturnArray.push_back( ValueArray[i] );
    }

    return ReturnArray;
}

void GetFirstRingVertices(
    ModelPart::MeshType& rInputMesh,
    std::map<IndexType, std::vector<IndexType>>& rOutputMap,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_verts,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_faces,
    OpenSubdiv::Far::TopologyRefiner* refiner = nullptr,
    IndexType ref_level = 0
    ) 
{
    /*
    Numbering based on Stam (example for regular vertex)

        9------2----->3
        ^      ^      |
        |      |      v
        8------1------4
        ^      |      |
        |      |      v
        7<-----6<-----5

    */
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices") << std::endl;
    std::map<IndexType, std::vector<IndexType>> first_ring_node_neighbourhoods;

    // if refiner not passed, generate a new one based on rInputMesh
    if (!refiner) {
        // std::map<IndexType, IndexType> map_verts, map_conds;
        refiner = GenerateOpenSubdivRefiner(rInputMesh, r_global_map_verts, r_global_map_faces, false);
        // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: no refiner specified... creating new one") << std::endl;
    }

    OpenSubdiv::Far::TopologyLevel const & refLevel = refiner->GetLevel(ref_level);
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: current refiner level ") << ref_level << std::endl;

    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: rInputMesh.NumberOfNodes() ") << rInputMesh.NumberOfNodes() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: refLevel.GetNumVertices() ") << refLevel.GetNumVertices() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: refLevel.GetNumFaces() ") << refLevel.GetNumFaces() << std::endl;
    SizeType num_vertices = refLevel.GetNumVertices();
    SizeType num_nodes = rInputMesh.NumberOfNodes();
    KRATOS_ERROR_IF_NOT(num_nodes == num_vertices);

    std::map<IndexType, IndexType> map_verts = r_global_map_verts[ref_level];
    std::map<IndexType, IndexType> map_faces = r_global_map_faces[ref_level];

    // initializes the 1-ring node neighbourhoods for every node in the input mesh
    
    // loop over all nodes to find neighbourhood of every node
    auto t1_v = high_resolution_clock::now();
    const auto nodes_begin = rInputMesh.NodesBegin();
#pragma omp parallel for
    for (IndexType vert_index = 0; vert_index < num_vertices; ++vert_index) {
        // find control node in mesh rInputMesh
        auto node_it = nodes_begin + vert_index;
        IndexType kratos_center_node_id = node_it->GetId();
        IndexType osd_center_vert_idx = GetOSDIndexFromKratosID(kratos_center_node_id, map_verts);
        // IndexType osd_center_vertex_index = r_center_node_id - 1;
        // skip free edge nodes since they need to be evaluated seperately
        if (refLevel.IsVertexBoundary(osd_center_vert_idx)) {
            #pragma omp critical
            rOutputMap[kratos_center_node_id] = {};
            continue;
        }

        std::vector<IndexType> temp_neighbourhood;
        temp_neighbourhood.push_back(kratos_center_node_id);
        ConstIndexArray face_indices = refLevel.GetVertexFaces(osd_center_vert_idx);
        
        // loop over every face that is adjacent to the vertex
        // (center_node and last index in temp_neighbourhood)
        IndexType starting_osd_index = osd_center_vert_idx;
        IndexType face_loop_count = 0;
        // for (SizeType j = 0; j < r_elements.size(); j++) {
        for (auto f_idx_it = face_indices.begin(); f_idx_it != face_indices.end(); ++f_idx_it) {

            IndexType f_idx = *f_idx_it;
            ConstIndexArray verts_indices = refLevel.GetFaceVertices(f_idx);

            // reorder the IDs such that the first ID is r_center_node_id
            std::vector<IndexType> reordered_verts_indices = ReorderIndicesToStartAtValue(starting_osd_index, verts_indices);
            // loop over the nodes of the element and identify if it has been added to the neighbourhood
            for(auto missing_it = reordered_verts_indices.begin()+1; missing_it != reordered_verts_indices.end(); ++missing_it) {
                
                IndexType kratos_node_idx = map_verts[*missing_it];
                bool is_present = (std::find(temp_neighbourhood.begin(), temp_neighbourhood.end(), kratos_node_idx) != temp_neighbourhood.end());

                if (is_present) {
                    break; // break out of loop, since necessary nodes from reordered_verts_indices are already present in neighbourhood
                }
                temp_neighbourhood.push_back(kratos_node_idx);
                starting_osd_index = *missing_it;
                // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: debug 5 : *missing_it ") << *missing_it << std::endl;
            }
        }
        // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: temp_neighbourhood.size() ") << temp_neighbourhood.size() << std::endl;
        #pragma omp critical
        rOutputMap[kratos_center_node_id] = temp_neighbourhood;
        // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: rOutputMap[r_center_node_id] ") << rOutputMap[r_center_node_id] << std::endl;
    }
    auto t2_v = high_resolution_clock::now();
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: number of vertices ") << num_vertices << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: time needed for loop ") << duration_cast<seconds>(t2_v-t1_v).count() << std::endl;

    // KRATOS_INFO("CATMULL_CLARK :: GetFirstRingVertices :: rOutputMap.size() ") << rOutputMap.size() << std::endl;
    // for(const auto& elem : rOutputMap)
    // {
    //     KRATOS_INFO("CATMULL_CLARK                 :: ") << elem.first << " : " << elem.second << std::endl;
    // }
};

bool GetSecondRingVertices(
    IndexType KratosFaceId,     // could be changed to OsdFaceIdx to avoid mapping two times for no reason
    std::vector<IndexType>& SecondRingVertices, // contains osd ids
    std::map<IndexType, std::map<IndexType, IndexType>>& r_map_verts,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_map_faces,
    OpenSubdiv::Far::TopologyRefiner* refiner = nullptr,
    IndexType refining_depth = 0
    )
{
    SecondRingVertices.clear();
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices") << std::endl;
    // assuming only quadrilaterals, since at least 1 subdivision has been performed previously
    bool is_irregular_patch = false;

    // std::map<IndexType, std::map<IndexType, IndexType>> map_vert_to_node, map_face_to_cond;
    // OpenSubdiv::Far::TopologyRefiner * refiner_for_geometry_description = GenerateOpenSubdivRefiner(rInputMesh, map_vert_to_node, map_face_to_cond, false);

    // get local patch geometry around osd_face_id
    // std::vector<int> num_vertices_per_face;
    // std::vector<int> vertex_indices_per_face;
    // int ref_level = refining_depth;
    // if (refining_depth == 0) ref_level = refiner->GetNumLevels() - 1;
    // else ref_level = refining_depth;

    int available_ref = refiner->GetNumLevels() - 1;
    OpenSubdiv::Far::TopologyLevel const & refLevel = refiner->GetLevel(refining_depth);
    KRATOS_ERROR_IF(refining_depth > available_ref) << "refining depth not available" << std::endl;
    // int ref_level = 0;
    // OpenSubdiv::Far::TopologyLevel const & refLevel = refiner_for_geometry_description->GetLevel(ref_level);

    // KRATOS_INFO("map_vert_to_node");
    // for (auto& it : map_vert_to_node[0]) {
    //     KRATOS_INFO(" ") << it.first << "->" << it.second << " # ";
    // }
    // KRATOS_INFO(" ") << std::endl;
    // KRATOS_INFO("map_face_to_cond");
    // for (auto& it : map_face_to_cond[0]) {
    //     KRATOS_INFO(" ") << it.first << "->" << it.second << " # ";
    // }
    // KRATOS_INFO(" ") << std::endl;
    // KRATOS_INFO("r_map_verts");
    // for (auto& it : r_map_verts[ref_level]) {
    //     KRATOS_INFO(" ") << it.first << "->" << it.second << " # ";
    // }
    // KRATOS_INFO(" ") << std::endl;
    // KRATOS_INFO("r_map_faces");
    // for (auto& it : r_map_faces[ref_level]) {
    //     KRATOS_INFO(" ") << it.first << "->" << it.second << " # ";
    // }
    // KRATOS_INFO(" ") << std::endl;

    // bool is_same = (map_vert_to_node[0] == r_map_verts[refining_depth]);
    // KRATOS_INFO("is_same") << is_same << std::endl;

    // KRATOS_INFO("refining_depth") << refining_depth << std::endl;
    // KRATOS_INFO("r_map_faces.size") << r_map_faces.size() << std::endl;
    // KRATOS_INFO("r_map_faces[refining_depth].size") << r_map_faces[refining_depth].size() << std::endl;

    IndexType osd_face_id = GetOSDIndexFromKratosID(KratosFaceId, r_map_faces[refining_depth]);
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: KratosFaceId") << KratosFaceId << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: osd_face_id") << osd_face_id << std::endl;

    // second ring vertices and faces
    ConstIndexArray vert_indices_of_face = refLevel.GetFaceVertices(osd_face_id);
    // vert_indices_of_face.size() should be 4 since there has been at least 1 subdivision previously

    // check if irregular
    for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
        // KRATOS_INFO("GetSecondRingVertices :: !refLevel.IsVertexValenceRegular(vert_indices_of_face[vert_i])") << !refLevel.IsVertexValenceRegular(vert_indices_of_face[vert_i]) << std::endl;
        if (!refLevel.IsVertexValenceRegular(vert_indices_of_face[vert_i])) is_irregular_patch = true;
        else if (refLevel.IsVertexBoundary(vert_indices_of_face[vert_i])) is_irregular_patch = true;
        // else if (refLevel.GetVertexFaces(vert_i).size() != 4) is_irregular_patch = true;
        
        ConstIndexArray faces_around_vert = refLevel.GetVertexFaces(vert_indices_of_face[vert_i]);
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: faces_around_vert.size") << faces_around_vert.size() << std::endl;
        for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                ConstIndexArray vertices_of_face_j = refLevel.GetFaceVertices(face_index);
                // KRATOS_INFO("GetSecondRingVertices :: vertices_of_face_j.size() != 4") << (vertices_of_face_j.size() != 4) << std::endl;
                if (vertices_of_face_j.size() != 4) is_irregular_patch = true;
        }
    }
    
    // create 2-ring neighbourhood
    if (is_irregular_patch) {
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: irregular patch ...") << std::endl;
        for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
            IndexType vert_index = vert_indices_of_face[vert_i];
            ConstIndexArray faces_around_vert = refLevel.GetVertexFaces(vert_index);

            for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                ConstIndexArray vertices_of_face_j = refLevel.GetFaceVertices(face_index);
                for (int vert_k = 0; vert_k < vertices_of_face_j.size(); ++vert_k) {
                    IndexType neighb_vert = vertices_of_face_j[vert_k];
                    auto v_it = std::find(SecondRingVertices.begin(), SecondRingVertices.end(), neighb_vert);
                    // add vertex index if not yet included
                    if (v_it == SecondRingVertices.end()) SecondRingVertices.push_back(neighb_vert);
                }
            }
        }
    }
    else {
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: regular patch ...") << std::endl;
        // prepare vector structure for regular patch
        SecondRingVertices.resize(16);

        std::vector<std::vector<IndexType>> patchPointsIndicesPerFace = { {  5,  4,  0,  1 },
                                                                        {  6,  2,  3,  7 },
                                                                        { 10, 11, 15, 14 },
                                                                        {  9, 13, 12,  8 } };
        
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: Found regular patch, creating ordered patch... ") << std::endl;
        // create ordered regular patch
        IndexType patch_position_index = 0;
        for (auto vert_index_it = vert_indices_of_face.begin(); vert_index_it != vert_indices_of_face.end(); ++vert_index_it) {
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: debug 1") << std::endl;
            IndexType vert_index = *vert_index_it;
            std::vector<IndexType> patch_indices = patchPointsIndicesPerFace[patch_position_index];
            patch_position_index += 1;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: debug 2") << std::endl;

            ConstIndexArray face_indices = refLevel.GetVertexFaces(vert_index);
            IndexType opposite_face_index = GetOppositeIndex(face_indices, osd_face_id);
            ConstIndexArray vert_indices_for_patch = refLevel.GetFaceVertices(opposite_face_index);
            std::vector<IndexType> reordered_vert_indices = ReorderIndicesToStartAtValue(vert_index, vert_indices_for_patch);

            // KRATOS_INFO("GetSecondRingVertices :: vert_indices_for_patch.size() == 4") << (vert_indices_for_patch.size() == 4) << std::endl;
            // KRATOS_INFO("GetSecondRingVertices :: reordered_vert_indices.size() == 4") << (reordered_vert_indices.size() == 4) << std::endl;
            for (int i = 0; i < 4; ++i) {
                SecondRingVertices[patch_indices[i]] = reordered_vert_indices[i];
            }
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: SecondRingVertices ") << SecondRingVertices << std::endl;
        }
    }
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: SecondRingVertices osd indices") << SecondRingVertices << std::endl;
    std::vector<IndexType> SecondRingVerticesKratosIds;
    for (int i = 0; i < SecondRingVertices.size(); ++i) {
        SecondRingVerticesKratosIds.push_back(r_map_verts[refining_depth][SecondRingVertices[i]]);
    }
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: SecondRingVertices kratos ids") << SecondRingVerticesKratosIds << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVertices :: Reached end of function") << std::endl;

    return is_irregular_patch;
}

bool GetSecondRingFaces(
    IndexType KratosFaceId, 
    std::vector<IndexType>& SecondRingFaces, // contains osd ids
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_verts,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_faces,
    OpenSubdiv::Far::TopologyRefiner* refiner = nullptr,
    IndexType refining_depth = 0
    )
{
    SecondRingFaces.clear();

    bool is_irregular_patch = false;

    int available_ref = refiner->GetNumLevels() - 1;
    OpenSubdiv::Far::TopologyLevel const & refLevel = refiner->GetLevel(refining_depth);
    KRATOS_ERROR_IF(refining_depth > available_ref) << "refining depth not available" << std::endl;
    
    IndexType osd_face_id = GetOSDIndexFromKratosID(KratosFaceId, r_global_map_faces[refining_depth]);

    // second ring faces
    SecondRingFaces.push_back(osd_face_id);
    ConstIndexArray vert_indices_of_face = refLevel.GetFaceVertices(osd_face_id);
    
    // check if irregular
    for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
        if (!refLevel.IsVertexValenceRegular(vert_indices_of_face[vert_i])) is_irregular_patch = true;
        
        ConstIndexArray faces_around_vert = refLevel.GetVertexFaces(vert_indices_of_face[vert_i]);
        for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                ConstIndexArray vertices_of_face_j = refLevel.GetFaceVertices(face_index);
                if (vertices_of_face_j.size() != 4) is_irregular_patch = true;
        }
    }

    // create 2-ring neighbourhood
    if (is_irregular_patch) {
        // KRATOS_INFO("GetSecondRingFaces :: irregular patch ...") << std::endl;

        for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
            IndexType vert_index = vert_indices_of_face[vert_i];
            ConstIndexArray faces_around_vert = refLevel.GetVertexFaces(vert_index);

            for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                // skip if face already accounted for, else add face to vector
                auto f_it = std::find(SecondRingFaces.begin(), SecondRingFaces.end(), face_index);
                if (f_it == SecondRingFaces.end()) SecondRingFaces.push_back(face_index);
            }
        }
    }
    else {
        // I guess else faces are not needed. 
        // They should only be needed to create the irregular patch in order to refine locally.
    }
    return is_irregular_patch;
}

std::vector<fcpw::Vector3i> TriangulateMesh(ModelPart::MeshType& rInputMesh, std::map<IndexType, IndexType>& MapTriaIndicesToConditions, std::map<IndexType, IndexType>& MapNodeIdsToIndices)
{
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: TriangulateMesh") << std::endl;

    std::vector<fcpw::Vector3i> triangle_indices;
    SizeType tria_index = 0;
    for (auto cond_it = rInputMesh.ConditionsBegin(); cond_it != rInputMesh.ConditionsEnd(); ++cond_it) {
        ConditionType& cond = *cond_it;
        GeometryType& geom = cond.GetGeometry();
        if (geom.size() == 3) {
            fcpw::Vector3i indices;
            indices << MapNodeIdsToIndices[ geom[0].GetId() ], MapNodeIdsToIndices[ geom[1].GetId() ], MapNodeIdsToIndices[ geom[2].GetId() ];
            triangle_indices.push_back(indices);
            MapTriaIndicesToConditions[tria_index] = cond.GetId();
            tria_index += 1;
        }
        else if (geom.size() == 4) {

            fcpw::Vector3i indices_1, indices_2;

            indices_1 << MapNodeIdsToIndices[ geom[0].GetId() ], MapNodeIdsToIndices[ geom[1].GetId() ], MapNodeIdsToIndices[ geom[2].GetId() ];
            indices_2 << MapNodeIdsToIndices[ geom[0].GetId() ], MapNodeIdsToIndices[ geom[2].GetId() ], MapNodeIdsToIndices[ geom[3].GetId() ];

            triangle_indices.push_back(indices_1);
            triangle_indices.push_back(indices_2);

            MapTriaIndicesToConditions[tria_index] = cond.GetId();
            MapTriaIndicesToConditions[tria_index + 1] = cond.GetId();

            tria_index += 2;
        }
        else if (geom.size() > 4) {
            KRATOS_ERROR << "Polygons with more than 4 vertices are not supported yet!" << std::endl;
        }
    }
    return triangle_indices;
}

IndexType FindClosestConditionOnMesh(NodeType& rFeNode, ModelPart::MeshType& rSearchMesh)
{
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh") << std::endl;

    // initialize a 3d scene
    fcpw::Scene<3> scene;

    std::vector<fcpw::Vector3> mesh_nodes;
    IndexType node_index = 0;
    std::map<IndexType, IndexType> map_node_id_to_index;
    // KRATOS_INFO(" -- ") << rSearchMesh.GetNode(1).Coordinates()[0] << rSearchMesh.GetNode(1).Coordinates()[1] << rSearchMesh.GetNode(1).Coordinates()[2] << std::endl;
    for (auto node_it = rSearchMesh.NodesBegin(); node_it != rSearchMesh.NodesEnd(); ++node_it) {
        NodeType& node = *node_it;
        fcpw::Vector3 new_node;
        new_node << node.Coordinates()[0], node.Coordinates()[1], node.Coordinates()[2];
        // KRATOS_INFO(" -- ") << node.Coordinates()[0] << ", " << node.Coordinates()[1] << ", " << node.Coordinates()[2] << ", " << std::endl;
        mesh_nodes.push_back(new_node);
        map_node_id_to_index[node.GetId()] = node_index;
        node_index += 1;
    }
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: before triangulation") << std::endl;

    std::map<IndexType, IndexType> map_tria_index_to_condition_id;
    const std::vector<fcpw::Vector3i>& triangle_indices = TriangulateMesh(rSearchMesh, map_tria_index_to_condition_id, map_node_id_to_index);
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: after triangulation") << std::endl;

    // load positions and indices of a single triangle mesh
    scene.setObjectCount(1);
    scene.setObjectVertices(mesh_nodes, 0);
    scene.setObjectTriangles(triangle_indices, 0);

    // build acceleration structure
    fcpw::AggregateType aggregateType = fcpw::AggregateType::Bvh_SurfaceArea;
    bool buildVectorizedBvh = true;
    scene.build(aggregateType, buildVectorizedBvh);
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: after fcpw build") << std::endl;

    // perform a closest point query
    fcpw::Interaction<3> interaction;
    fcpw::Vector3 queryPoint;
    queryPoint << rFeNode.Coordinates()[0], rFeNode.Coordinates()[1], rFeNode.Coordinates()[2];
    bool found = scene.findClosestPoint(queryPoint, interaction);
    // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: after closest point query") << std::endl;

    // access distance and closest point via interaction.d and interaction.p (resp.)
    if (found) {
        int fcpw_tria_index = interaction.primitiveIndex;
        // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: map_tria_index_to_condition_id.size ") << map_tria_index_to_condition_id.size() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: FindClosestConditionOnMesh :: Closest Condition is KratosFaceId ") << map_tria_index_to_condition_id[fcpw_tria_index] << std::endl;
        return map_tria_index_to_condition_id[fcpw_tria_index];
    }
    else KRATOS_ERROR << "FindClosestConditionOnMesh :: Could not determine nearest point on tria surface!" << std::endl;
}

std::vector<double> GetLimitWeights(SizeType NeighbourhoodSize)
{
    // The weights for computing the limit stem from an analytical approach 
    // where the spectral analysis of the local subdivision matrix is used to 
    // determine the weights in the limit of every point in the neighbourhood.
    // This enables to compute every point on the limit, even irregular points.
    if (NeighbourhoodSize == 7)     // valence 3
        {
            return {0.375, 1.0/6.0, 1.0/24.0, 1.0/6.0, 1.0/24.0, 1.0/6.0, 1.0/24.0};
        }
        else if (NeighbourhoodSize == 9)    // valence 4, regular point for Catmull-Clark
        {
            return {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};
        }
        else if (NeighbourhoodSize == 11)   //valence 5
        {
            return {0.5, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02};
        }
        else if (NeighbourhoodSize == 13)   // valence 6
        {
            return {6.0/11.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 
            5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0};
        }
        else if (NeighbourhoodSize == 15)   // valence 7
        {
            return {7.0/12.0, 0.04761905, 0.01190476, 0.04761905, 0.01190476, 0.04761905, 0.01190476, 
            0.04761905, 0.01190476, 0.04761905, 0.01190476, 0.04761905, 0.01190476, 0.04761905, 0.01190476};
        }
        else if (NeighbourhoodSize == 17)   // valence 8
        {
            return {0.61538462, 0.03846154, 0.00961538, 0.03846154, 0.00961538, 0.03846154, 0.00961538, 0.03846154, 0.00961538, 
            0.03846154, 0.00961538, 0.03846154, 0.00961538, 0.03846154, 0.00961538, 0.03846154, 0.00961538};
        }
        else if (NeighbourhoodSize == 19)   // valence 9
        {
            return {0.64285714, 0.03174603, 0.00793651, 0.03174603, 0.00793651, 0.03174603, 0.00793651, 0.03174603, 0.00793651, 
            0.03174603, 0.00793651, 0.03174603, 0.00793651, 0.03174603, 0.00793651, 0.03174603, 0.00793651, 0.03174603, 0.00793651};
        }        
        else { // NeighbourhoodSize > 19 not implemented
            // TODO: eventuallly fall back to computing limit map analytically
            KRATOS_ERROR << "Neighbourhood size for pushing nodes to limit is too large (max. is 19 with valence 9). Limit map is not implemented for neighbourhood size " << NeighbourhoodSize << std::endl;
        }
}

void PushNodesToLimit(ModelPart& rInputModelPart, ModelPart& rOutputModelPart, std::map<IndexType, std::vector<IndexType>> rFirstRingNodeNeighbourhoods)
{   
    auto t1_v = high_resolution_clock::now();
    // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit") << std::endl;
    const auto nodes_begin = rInputModelPart.NodesBegin();
// #pragma omp parallel for
    for (IndexType node_idx = 0; node_idx < rInputModelPart.NumberOfNodes(); ++node_idx)
    {
        auto node_it = nodes_begin + node_idx;
        IndexType node_id = node_it->GetId();
        // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit :: Current Kratos node ID") << node_id << std::endl;
        // const CoordinatesArrayType& r_coords = r_node_i.Coordinates();

        std::vector<IndexType> r_node_i_neighbour_node_ids = rFirstRingNodeNeighbourhoods[node_id];
        // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit :: with first ring neighbourhood") << r_node_i_neighbour_node_ids << std::endl;
        
        SizeType neighb_size = r_node_i_neighbour_node_ids.size();
        // skip boundary nodes as they cannot be mapped to limit
        // boundary nodes do not have a 1-ring neighb (it was not created)
        if (neighb_size == 0) {
            // remove node from limit model part
            rOutputModelPart.RemoveNode(rInputModelPart.pGetNode(node_id));
            continue;
        }

        std::vector<CoordinatesArrayType> neighbourhood_coords;
        
        for (auto neighb_idx_it = r_node_i_neighbour_node_ids.begin(); neighb_idx_it != r_node_i_neighbour_node_ids.end(); ++neighb_idx_it)
        {
            neighbourhood_coords.push_back(rInputModelPart.GetNode(*neighb_idx_it).Coordinates());
        }

        std::vector<double> limit_weights = GetLimitWeights(neighb_size);
        std::vector<double> coords = {0.0, 0.0, 0.0};
        for(SizeType k = 0; k < limit_weights.size(); k++) {
            double weight = limit_weights[k];
            coords[0] += weight * neighbourhood_coords[k].operator()(0);
            coords[1] += weight * neighbourhood_coords[k].operator()(1);
            coords[2] += weight * neighbourhood_coords[k].operator()(2);
        }
        rOutputModelPart.GetNode(node_id).Coordinates()[0] = coords[0];
        rOutputModelPart.GetNode(node_id).Coordinates()[1] = coords[1];
        rOutputModelPart.GetNode(node_id).Coordinates()[2] = coords[2];
    }
    // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit successful") << std::endl;
    auto t2_v = high_resolution_clock::now();
    // KRATOS_INFO("OPENSUBDIV_UTILS :: PushNodesToLimit :: number of vertices ") << rOutputModelPart.NumberOfNodes() << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: PushNodesToLimit :: time needed for loop ") << duration_cast<seconds>(t2_v-t1_v).count() << std::endl;
}

std::vector<double> GetUnitWeightDistributionForNode(
    IndexType osd_v_idx, 
    OpenSubdiv::Far::TopologyRefiner * refiner, 
    SizeType ref_level
    )
{
    // Create a buffer to hold the position of the refined verts and
    // local points, then copy the coarse positions at the beginning.
    
    // SizeType n_uniformRef = refiner->GetNumLevels() - 1;
    SizeType nverts = refiner->GetLevel(0).GetNumVertices();
    // KRATOS_INFO("CATMULL_CLARK :: GetUnitWeightDistributionForNode: ref_level") << ref_level << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: GetUnitWeightDistributionForNode: nverts") << nverts << std::endl;

    // std::vector<double> vertices(3*nverts, 0.0);
    std::vector<double> cp_weights(nverts, 0.0);
    cp_weights[osd_v_idx] = 1.0;

    // Vertex *verts = (Vertex*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(Vertex) );
    // memcpy(&verts[0], &vertices[0], nverts*sizeof(double));

    VertexValue *vverts = (VertexValue*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(VertexValue) );
    memcpy(&vverts[0], &cp_weights[0], nverts*sizeof(double));

    // Interpolate vertex primvar data at each level
    OpenSubdiv::Far::PrimvarRefiner primvarRefiner(*refiner);

    // Vertex * src = verts;
    VertexValue * vsrc = vverts;
    for (SizeType level = 1; level <= ref_level; ++level) {
        // Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
        // primvarRefiner.Interpolate(level, src, dst);
        // src = dst;

        VertexValue * vdst = vsrc + refiner->GetLevel(level-1).GetNumVertices();
        primvarRefiner.Interpolate(level, vsrc, vdst);
        vsrc = vdst;
    }

    // return the current vertices and faces
    OpenSubdiv::Far::TopologyLevel const & refLastLevel = refiner->GetLevel(ref_level);
    
    nverts = refLastLevel.GetNumVertices();
    
    // copy vertex information
    int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
    // vertices.resize(nverts*3);
    // memcpy(&vertices[0], &verts[firstOfLastVerts], nverts*3*sizeof(double));

    cp_weights.resize(nverts);
    memcpy(&cp_weights[0], &vverts[firstOfLastVerts], nverts*sizeof(double));

    // KRATOS_INFO("CATMULL_CLARK :: GetUnitWeightDistributionForNode :: cp_weights computed") << cp_weights << std::endl; // << cp_weights << std::endl;
    // refLastLevel.GetNumVertices
    // free(verts);
    // free(vverts);
    std::vector<double> cp_weights_copy = cp_weights;   // should be deep copy
    return cp_weights_copy;
}

void OutputModelPartToVtk(ModelPart& rOuptutMP, std::string filename)
{
    // output for debug only
    Parameters parameters = Parameters(R"(
    {    
        "model_part_name"                             : "refined_patch",
        "file_format"                                 : "ascii",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_deformed_configuration"                : true,
        "write_ids"                                   : false,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : [],
        "entity_type"                                 : "automatic",
        "output_path"                                 : "VTK_Output"
    })" );

    VtkOutput(rOuptutMP, parameters).PrintOutput(filename);
}

void FillPatchNodesVector(
    std::vector<NodeTypePointer>& rPatchNodesVector, 
    ModelPart::MeshType& rMeshGrid,
    std::vector<IndexType>& rSecondRingVertices,
    std::map<IndexType, IndexType>& rMapVertToNode
    )
{
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FillPatchNodeVector") << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rSecondRingVertices") << rSecondRingVertices << std::endl;
    for (IndexType i = 0; i < rSecondRingVertices.size(); ++i) {
        IndexType kratos_id = rMapVertToNode[rSecondRingVertices[i]];
        // KRATOS_INFO("rSecondRingVertices[i]") << rSecondRingVertices[i] << " ... kratos_id: " << kratos_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FillPatchNodeVector :: rSecondRingVertices[i] = osd_id ") << rSecondRingVertices[i] << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FillPatchNodeVector :: kratos_id ") << kratos_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FillPatchNodeVector :: rMeshGrid.GetNode(kratos_id) ") << rMeshGrid.GetNode(kratos_id) << std::endl;
        rPatchNodesVector.push_back(rMeshGrid.pGetNode(kratos_id));
    }
}

void GetRefiningWeights(
    std::map< IndexType, std::map<IndexType, double> >& rRefiningWeights,
    std::vector<IndexType>& rBaseSecondRingVertices,
    std::vector<IndexType>& rPatchIds,
    std::map< IndexType, std::map<IndexType, IndexType> >& rGlobalMapVertToNode,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr,
    SizeType refining_depth = 1
)
{
    for(auto v_it = rBaseSecondRingVertices.begin(); v_it != rBaseSecondRingVertices.end(); ++v_it) {
        std::map<IndexType, double> influence_weights_of_base_on_refined;   // size is 16 if patch is regular
        IndexType osd_base_v_id = *v_it;
        IndexType kratos_base_v_id = rGlobalMapVertToNode[0][osd_base_v_id];
        
        std::vector<double> weight_distr_of_cp_i_on_refined_patch = GetUnitWeightDistributionForNode(osd_base_v_id, refiner, refining_depth);


        // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: rPatchIds ") << rPatchIds << std::endl;
        IndexType loop_nr = 0;
        for(auto regular_patch_it = rPatchIds.begin(); regular_patch_it != rPatchIds.end(); ++ regular_patch_it) {
            IndexType osd_patch_vertex_id = *regular_patch_it;
            IndexType kratos_patch_node_id = rGlobalMapVertToNode[refining_depth][osd_patch_vertex_id];
            influence_weights_of_base_on_refined[kratos_patch_node_id] = weight_distr_of_cp_i_on_refined_patch[osd_patch_vertex_id];
            // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: loop_nr ") << loop_nr << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: osd_patch_vertex_id ") << osd_patch_vertex_id << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: kratos_patch_node_id ") << kratos_patch_node_id << std::endl;
            ++loop_nr;
        }
        // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: refining_depth ") << refining_depth << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: kratos_base_v_id ") << kratos_base_v_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetRefiningWeights :: GetRefiningWeights :: influence_weights_of_base_on_refined ") << influence_weights_of_base_on_refined << std::endl;

        rRefiningWeights[kratos_base_v_id] = influence_weights_of_base_on_refined;
    }
}

void SetUpPatchMP(
    ModelPart& rGrid,
    ModelPart& rPatchMP, 
    std::vector<IndexType>& rSecondRingVertices,
    std::vector<IndexType>& rSecondRingFaces,
    std::map<IndexType, IndexType>& r_map_vert_to_node,
    std::map<IndexType, IndexType>& r_map_face_to_cond
)
{
    // KRATOS_INFO("CATMULL_CLARK :: SetUpPatchMP") << std::endl;

    SizeType nverts = rSecondRingVertices.size();
    NodesContainerType::Pointer p_new_nodes_container = make_shared<NodesContainerType>();
    // p_new_nodes_container->reserve(nverts);
    for (IndexType i = 0; i < nverts; ++i) {
        IndexType osd_vert_idx = rSecondRingVertices[i];
        IndexType kratos_node_id = r_map_vert_to_node[osd_vert_idx];
        CoordinatesArrayType node_coords = rGrid.GetNode(kratos_node_id).Coordinates();
        std::vector<double> new_node_coord(3);
        new_node_coord[0] = node_coords[0];
        new_node_coord[1] = node_coords[1];
        new_node_coord[2] = node_coords[2];
        NodeTypePointer p_node = make_intrusive<Node>(kratos_node_id, new_node_coord);
        p_new_nodes_container->push_back(p_node);
    }

    p_new_nodes_container->Unique();

    rPatchMP.SetNodes(p_new_nodes_container);
    rPatchMP.GetCommunicator().LocalMesh().SetNodes(p_new_nodes_container);

    SizeType nfaces = rSecondRingFaces.size();
    std::string condition_name = "SurfaceCondition3D4N";
    
    Properties::Pointer p_property;
    if (rPatchMP.HasProperties(0)) {
        p_property = rPatchMP.pGetProperties(0);
    }
    else {
        p_property = rPatchMP.CreateNewProperties(0);
    }
    for (IndexType i = 0; i < nfaces; ++i) {
        IndexType osd_face_idx = rSecondRingFaces[i];
        IndexType kratos_cond_id = r_map_face_to_cond[osd_face_idx];
        const ConditionType& r_cond = rGrid.GetCondition(kratos_cond_id);
        const GeometryType& r_geom = r_cond.GetGeometry();
        SizeType number_of_nodes = r_geom.PointsNumber();
        std::vector<IndexType> node_id_array(number_of_nodes);
        for (IndexType j = 0; j < number_of_nodes; ++j) {
            node_id_array[j] = r_geom[j].GetId();
        }
        rPatchMP.CreateNewCondition(condition_name, kratos_cond_id, node_id_array, p_property);
    }

    rPatchMP.GetCommunicator().LocalMesh().SetConditions(rPatchMP.pConditions());
    
}

void CopyModelPart(const ModelPart& rCopyMP, ModelPart& rDestinationMP) 
{
    // KRATOS_INFO("CATMULL_CLARK :: CopyModelPart") << std::endl;


    SizeType nverts = rCopyMP.NumberOfNodes();
    const auto& nodes_begin = rCopyMP.NodesBegin();
    NodesContainerType::Pointer p_new_nodes_container = make_shared<NodesContainerType>();
    // p_new_nodes_container->reserve(nverts);
    for (IndexType i = 0; i < nverts; ++i) {
        const auto& node_it = nodes_begin + i;
        IndexType kratos_node_id = node_it->GetId();
        CoordinatesArrayType node_coords = node_it->Coordinates();
        std::vector<double> new_node_coord(3);
        new_node_coord[0] = node_coords[0];
        new_node_coord[1] = node_coords[1];
        new_node_coord[2] = node_coords[2];
        NodeTypePointer p_node = make_intrusive<Node>(kratos_node_id, new_node_coord);
        p_new_nodes_container->push_back(p_node);
    }

    p_new_nodes_container->Unique();

    rDestinationMP.SetNodes(p_new_nodes_container);
    rDestinationMP.GetCommunicator().LocalMesh().SetNodes(p_new_nodes_container);

    SizeType nconds = rCopyMP.NumberOfConditions();
    const auto& conds_begin = rCopyMP.ConditionsBegin();
    Properties::Pointer p_property;
    if (rDestinationMP.HasProperties(0)) {
        p_property = rDestinationMP.pGetProperties(0);
    }
    else {
        p_property = rDestinationMP.CreateNewProperties(0);
    }
    for (IndexType i = 0; i < nconds; ++i) {
        const auto& cond_it = conds_begin + i;
        IndexType kratos_cond_id = cond_it->GetId();
        const GeometryType& r_geom = cond_it->GetGeometry();
        SizeType number_of_nodes = r_geom.PointsNumber();

        std::string condition_name = "SurfaceCondition3D4N";
        if (number_of_nodes == 3) {
            condition_name = "SurfaceCondition3D3N";
        }
        else if (number_of_nodes > 4) {
            KRATOS_ERROR << "Polygons with more than 4 nodes are not supported yet!" << std::endl;
        }

        std::vector<IndexType> node_id_array(number_of_nodes);
        for (IndexType j = 0; j < number_of_nodes; ++j) {
            node_id_array[j] = r_geom[j].GetId();
        }
        
        // KRATOS_INFO("number_of_nodes") << number_of_nodes << std::endl;
        // KRATOS_INFO("condition_name") << condition_name << std::endl;
        // KRATOS_INFO("kratos_cond_id") << kratos_cond_id << std::endl;
        // KRATOS_INFO("node_id_array") << node_id_array << std::endl;
        // KRATOS_INFO("p_property") << p_property << std::endl;
        
        rDestinationMP.CreateNewCondition(condition_name, kratos_cond_id, node_id_array, p_property);
    }

    rDestinationMP.GetCommunicator().LocalMesh().SetConditions(rDestinationMP.pConditions());
    
}

std::vector<NodeTypePointer> CreateRegularPatchForFeNode(
    NodeType& rFeNode, 
    ModelPart& rControlPolygon,
    ModelPart& rGrid,
    ModelPart& rGridLimit,
    std::map<IndexType, Vector>& RefiningWeights,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_vert_to_node,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_face_to_cond,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode") << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rGrid.GetNode(1).Coordinates") << std::endl;
    // KRATOS_INFO(" -- ") << rGrid.GetNode(1).Coordinates()[0] << rGrid.GetNode(1).Coordinates()[1] << rGrid.GetNode(1).Coordinates()[2] << std::endl;
    SizeType max_refining_depth = 5;
    SizeType current_refining_depth = 0;    // we start at first refinement of control polygon
    Model patch_model;
    Model refined_patch_model;
    Model refined_limit_model;
    ModelPart& r_patch_mp = patch_model.CreateModelPart("patch", 3);
    ModelPart& r_refined_patch_mp = refined_patch_model.CreateModelPart("refined_patch", 3);
    ModelPart& r_refined_limit_mp = refined_limit_model.CreateModelPart("refined_limit", 3);
    r_patch_mp.CreateNewProperties(0);
    r_refined_patch_mp.CreateNewProperties(0);
    r_refined_limit_mp.CreateNewProperties(0);
    std::map<IndexType, std::vector<IndexType>> first_ring_vert_map;

    // 5.1 find closest control face on limit surface
    IndexType KratosFaceId = FindClosestConditionOnMesh(rFeNode, rGridLimit.GetMesh());
    IndexType FirstKratosFaceId = KratosFaceId;
    std::vector<IndexType> base_refined_second_ring_vertices;    // contains osd ids
    std::vector<IndexType> base_refined_second_ring_faces;       // contains osd ids
        
    // 5.2 create second ring vertices around control face
    bool patch_is_irregular = GetSecondRingVertices(
        KratosFaceId, 
        base_refined_second_ring_vertices, 
        r_global_map_vert_to_node,
        r_global_map_face_to_cond,
        refiner,
        1
    );
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refined_second_ring_vertices ") << base_refined_second_ring_vertices << std::endl;

    GetSecondRingFaces(
        KratosFaceId, 
        base_refined_second_ring_faces, 
        r_global_map_vert_to_node,
        r_global_map_face_to_cond,
        refiner,
        1
    );
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refined_second_ring_faces ") << base_refined_second_ring_faces << std::endl;

    // set up local patch to be refined
    SetUpPatchMP(
        rGrid, 
        r_patch_mp, 
        base_refined_second_ring_vertices, 
        base_refined_second_ring_faces, 
        r_global_map_vert_to_node[1], 
        r_global_map_face_to_cond[1]
        );
    // r_patch_mp.SetNodes(rInputMesh.pNodes());
    // r_patch_mp.SetConditions(rInputMesh.pConditions());
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;
    // ModelPart& r_patch_mp = AuxiliarModelPartUtilities(rInputModelPart).DeepCopyModelPart("patch", &patch_model);
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: after DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;

    // NodeTypePointer as return type
    std::vector<NodeTypePointer> rPatchNodesVector;

    // in case first refinement was enough, fill rPatchNodesVector and RefiningWeights, then return
    // still need to adapt refinement weights for the case where only one refinement is needed
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_is_irregular") << patch_is_irregular << std::endl;

    // optional: write files with refined patches
    std::string filename = r_patch_mp.Name() + "_ref" + std::to_string(current_refining_depth) + "_id" + std::to_string(rFeNode.GetId());
    r_patch_mp.Check();
    // OutputModelPartToVtk(r_patch_mp, filename);

    if (!patch_is_irregular) {

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: fill patch nodes vector with vertices from first refinement ") << std::endl;
        FillPatchNodesVector(rPatchNodesVector, rGrid.GetMesh(), base_refined_second_ring_vertices, r_global_map_vert_to_node[1]);

        // get weights from refining: meaning the influence of the base control points on the refined vertices of the regular patch
        IndexType osd_face_index = GetOSDIndexFromKratosID(KratosFaceId, r_global_map_face_to_cond[1]);
        OpenSubdiv::Far::TopologyLevel const & lastRefLevel = refiner->GetLevel(1);
        // get face index of control grid where FE node lies in
        for (SizeType ref = 0; ref < 1; ++ref) {
            osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
        }

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent osd_face_index on control polygon ") << osd_face_index << std::endl;
        IndexType kratos_base_face_index = r_global_map_face_to_cond[0][osd_face_index];
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent kratos_base_face_index on control polygon ") << kratos_base_face_index << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FirstKratosFaceId  ") << FirstKratosFaceId << std::endl;
        

        std::vector<IndexType> base_second_ring_vertices;    // contains osd ids
        // get second ring verts and faces for the base control grid P0
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: computing base_second_ring_vertices for refining weights due to 1 refinement ") << std::endl;
        GetSecondRingVertices(
            kratos_base_face_index, 
            base_second_ring_vertices, 
            r_global_map_vert_to_node, 
            r_global_map_face_to_cond,
            refiner,
            0
        );
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_second_ring_vertices ") << base_second_ring_vertices << std::endl;

        // now base_second_ring_vertices contains all osd ids of the ctrl pts that influence the vertices of final regular patch
        // evaluate unit weight distribution for every node on the base level (refinement 0)
        std::map<IndexType, std::map<IndexType, double> > refining_weights;
        GetRefiningWeights(refining_weights, base_second_ring_vertices, base_refined_second_ring_vertices, r_global_map_vert_to_node, refiner, 1);
        
        /// bring weights together
        const auto base_verts_begin = base_second_ring_vertices.begin();
        const auto base_refined_verts_begin = base_refined_second_ring_vertices.begin();
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: combine weights ") << std::endl;
        for (IndexType i = 0; i < base_second_ring_vertices.size(); ++i) {
            Vector weight_mult(16);
            weight_mult *= 0.0;
            auto base_v_osd_idx_it = base_verts_begin + i;
            IndexType base_kratos_node_id = r_global_map_vert_to_node[0][*base_v_osd_idx_it];

            for (IndexType k = 0; k < weight_mult.size(); ++k) {
                auto base_refined_verts_osd_idx_it = base_refined_verts_begin + k;
                IndexType regular_patch_kratos_node_id = r_global_map_vert_to_node[1][*base_refined_verts_osd_idx_it];
                weight_mult[k] += refining_weights[base_kratos_node_id][regular_patch_kratos_node_id];
            }

            RefiningWeights[base_kratos_node_id] = weight_mult;
        }

        return rPatchNodesVector;
    }
    // else if we have irregular patch:


    // refine the patch until fe-node is inside regular patch
    std::map<IndexType, std::map<IndexType, IndexType>> r_patch_map_vert_to_node;
    std::map<IndexType, std::map<IndexType, IndexType>> r_patch_map_face_to_cond;
    OpenSubdiv::Far::TopologyRefiner* patch_refiner = GenerateOpenSubdivRefiner(r_patch_mp.GetMesh(), r_patch_map_vert_to_node, r_patch_map_face_to_cond, false);

    std::vector<IndexType> regular_second_ring_vertices;    // regular patch
    while (patch_is_irregular && max_refining_depth > current_refining_depth) 
    {
        current_refining_depth += 1;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: current_refining_depth of patch") << current_refining_depth << std::endl;

        // auto geom = rInputMesh.GetCondition(KratosFaceId).GetGeometry();
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom ") << geom << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom IDs ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;
        RefineGeometry(r_patch_mp.GetMesh(), r_refined_patch_mp, r_patch_map_vert_to_node, r_patch_map_face_to_cond, current_refining_depth, patch_refiner);
    
        std::map<IndexType, std::vector<IndexType>> first_ring_vert_map;
        GetFirstRingVertices(r_refined_patch_mp.GetMesh(), first_ring_vert_map, r_patch_map_vert_to_node, r_patch_map_face_to_cond, patch_refiner, current_refining_depth);
        // create modelpart whose nodes lie on the limit surface
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << r_refined_patch_mp.Name() << std::endl;
        ModelPart& r_refined_limit_mp = AuxiliarModelPartUtilities(r_refined_patch_mp).DeepCopyModelPart(r_refined_patch_mp.Name()+"_limit_"+std::to_string(current_refining_depth), &refined_limit_model);
        PushNodesToLimit(r_refined_patch_mp, r_refined_limit_mp, first_ring_vert_map);
        
        // after refining: check which face fe_node lies in
        KratosFaceId = FindClosestConditionOnMesh(rFeNode, r_refined_limit_mp.GetMesh());
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New KratosFaceId in refined mesh ") << KratosFaceId << std::endl;

        // create 2-ring vertex neighbourhoods around the face the fe-node lies in
        patch_is_irregular = GetSecondRingVertices(
            KratosFaceId, 
            regular_second_ring_vertices, 
            r_patch_map_vert_to_node,
            r_patch_map_face_to_cond,
            patch_refiner,
            current_refining_depth
            );
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New patch for KratosFaceId: patch_is_irregular ") << patch_is_irregular << std::endl;
        // optional: write files with refined patches
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: r_refined_patch_mp.NumberOfNodes() ") << r_refined_patch_mp.NumberOfNodes() << std::endl;
        std::string filename = r_refined_patch_mp.Name() + "_ref" + std::to_string(current_refining_depth) + "_id" + std::to_string(rFeNode.GetId());
        r_refined_patch_mp.Check();
        // OutputModelPartToVtk(r_refined_patch_mp, filename);
    }
    // // optional: write files with refined patches
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: r_refined_patch_mp.NumberOfNodes() ") << r_refined_patch_mp.NumberOfNodes() << std::endl;
    // std::string filename = r_refined_patch_mp.Name() + "_ref" + std::to_string(current_refining_depth) + "_id" + std::to_string(rFeNode.GetId());
    // r_refined_patch_mp.Check();
    // OutputModelPartToVtk(r_refined_patch_mp, filename);

    if (current_refining_depth >= max_refining_depth && patch_is_irregular) {
        // KRATOS_INFO("Maximum refining depth reached, falling back to checking if fe-node is next to an irregular vertex ...") << std::endl;

        // actually near enoug to an irregular vertex to just take limit map of the neighbourhood
        // 1. get nearest node on base cp to fe-node (should be irregular vertex)
        OpenSubdiv::Far::TopologyLevel const & level = refiner->GetLevel(1);
        OpenSubdiv::Far::TopologyLevel const & firstLevel = refiner->GetLevel(0);
        IndexType OsdFaceIndexLevel1 = GetOSDIndexFromKratosID(FirstKratosFaceId, r_global_map_vert_to_node[1]);
        IndexType parent_face_osd_idx = level.GetFaceParentFace(OsdFaceIndexLevel1);
        // condition on level 1 that has irregular vertex:
        IndexType parent_cond_kratos_id = r_global_map_vert_to_node[0][parent_face_osd_idx];
        ConditionType original_cond = rControlPolygon.GetCondition(parent_cond_kratos_id);
        // get first ring vertices to push control polygon to limit
        std::map<IndexType, std::vector<IndexType>> base_first_ring_vert_map;   // contains only kratos node ids
        GetFirstRingVertices(rControlPolygon.GetMesh(), base_first_ring_vert_map, r_global_map_vert_to_node, r_global_map_face_to_cond, refiner, 0);
        // ModelPart& r_cp_limit_mp = AuxiliarModelPartUtilities(rControlPolygon).DeepCopyModelPart(rControlPolygon.Name()+"_limit");
        // Model& r_base_model = rControlPolygon.GetModel();
        Model limit_base_model;
        ModelPart& r_cp_limit_mp = limit_base_model.CreateModelPart("limit_control_points", 3);
        // r_cp_limit_mp.CreateNewProperties(0);
        // Mesh mesh_cp_clone = rControlPolygon.GetMesh().Clone();
        // r_cp_limit_mp.SetNodes(mesh_cp_clone.pNodes());
        // r_cp_limit_mp.SetConditions(mesh_cp_clone.pConditions());
        CopyModelPart(rControlPolygon, r_cp_limit_mp);
        // KRATOS_INFO("rControlPolygon.NumberOfNodes") << rControlPolygon.NumberOfNodes() << std::endl;
        // KRATOS_INFO("rControlPolygon.NumberOfConditions") << rControlPolygon.NumberOfConditions() << std::endl;
        // KRATOS_INFO("rControlPolygon.GetNode(1)") << rControlPolygon.GetNode(1) << std::endl;

        PushNodesToLimit(rControlPolygon, r_cp_limit_mp, base_first_ring_vert_map);
        // KRATOS_INFO("r_cp_limit_mp.NumberOfNodes") << r_cp_limit_mp.NumberOfNodes() << std::endl;
        // KRATOS_INFO("r_cp_limit_mp.NumberOfConditions") << r_cp_limit_mp.NumberOfConditions() << std::endl;
        // KRATOS_INFO("rControlPolygon.GetNode(1)") << rControlPolygon.GetNode(1) << std::endl;
        ConditionType original_cond_limit = r_cp_limit_mp.GetCondition(parent_cond_kratos_id);
        

        ConditionType& r_base_limit_kratos_cond = r_cp_limit_mp.GetCondition(parent_cond_kratos_id);
        GeometryType& r_geom = r_base_limit_kratos_cond.GetGeometry();
        CoordinatesArrayType& r_fe_coords = rFeNode.Coordinates();
        double min_dist = 1e12;
        IndexType closest_node_id;
        for (int i = 0; i < r_geom.PointsNumber(); ++i) {
            CoordinatesArrayType& cp_coords = r_geom[i].Coordinates();
            double dist = norm_2(r_fe_coords - cp_coords);
            if (dist < min_dist) {
                min_dist = dist;
                closest_node_id = r_geom[i].GetId();
            }
        }
        // check if closest_node_id is irregular node
        IndexType closest_node_osd_index = GetOSDIndexFromKratosID(closest_node_id, r_global_map_vert_to_node[0]);
        bool is_regular_valence = firstLevel.IsVertexValenceRegular(closest_node_osd_index);
        if (is_regular_valence) {
            // throw error
            KRATOS_ERROR << "Reached maximum refining depth, could not determine regular patch for FE-node " << rFeNode.GetId() << std::endl;
        }
        // else is irregular
        // KRATOS_INFO("FE-node distance to irregular vertex") << min_dist << std::endl;
        // 2. evaluate limit map in 1-ring neighb of irregular vertex to get weights of surrounding vertices
        std::vector<IndexType> first_ring_vertices_of_closest_node = base_first_ring_vert_map[closest_node_id];
        SizeType neighb_size = first_ring_vertices_of_closest_node.size();
        // KRATOS_INFO("neighb_size of irregular vertex") << neighb_size << std::endl;
        // get limit weights for irregular point
        std::vector<double> limit_weights = GetLimitWeights(neighb_size);
        Vector weights(neighb_size);
        for (IndexType i = 0; i < neighb_size; ++i) {
            weights[i] = limit_weights[i];
        }
        // KRATOS_INFO("Irregular Patch weights") << weights << std::endl;
        RefiningWeights[0] = weights;
        // FillPatchNodesVector(rPatchNodesVector, rControlPolygon.GetMesh(), first_ring_vertices_of_closest_node, r_global_map_vert_to_node[0]);
        for (IndexType i = 0; i < first_ring_vertices_of_closest_node.size(); ++i) {
            IndexType kratos_id = first_ring_vertices_of_closest_node[i];
            rPatchNodesVector.push_back(rControlPolygon.pGetNode(kratos_id));
        }
        // KRATOS_INFO("Found limit weights for irregular point!") << std::endl;
        // KRATOS_INFO("rPatchNodesVector") << rPatchNodesVector << std::endl;
        return rPatchNodesVector;
    }
    // else regular patch was found after refining 
    

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: fill patch nodes vector with vertices from refinement level ") << current_refining_depth << std::endl;
    FillPatchNodesVector(rPatchNodesVector, r_refined_patch_mp.GetMesh(), regular_second_ring_vertices, r_patch_map_vert_to_node[current_refining_depth]);

    // get weights from refining: meaning the influence of the base control points on the refined vertices of the regular patch
    /// 1. get weights for base level of patch_refiner
    IndexType osd_face_index = GetOSDIndexFromKratosID(KratosFaceId, r_patch_map_face_to_cond[current_refining_depth]);
    OpenSubdiv::Far::TopologyLevel const & lastRefLevel = patch_refiner->GetLevel(current_refining_depth);
    // get face index of control grid where FE node lies in
    for (SizeType ref = 0; ref < current_refining_depth; ++ref) {
        osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
    }

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent osd_face_index on control polygon ") << osd_face_index << std::endl;
    IndexType kratos_patch_base_face_index = r_patch_map_face_to_cond[0][osd_face_index];

    std::vector<IndexType> patch_base_second_ring_verts;    // contains osd ids
    // get second ring verts and faces for the base control grid P0
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: compute patch_base_second_ring_verts to get refining weights for refinement number ") << current_refining_depth << std::endl;
    GetSecondRingVertices(
        FirstKratosFaceId, 
        patch_base_second_ring_verts, 
        r_patch_map_vert_to_node, 
        r_patch_map_face_to_cond,
        patch_refiner,
        0
    );
    
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: loop over patch_base_second_ring_verts ") << patch_base_second_ring_verts << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_ring_kratos_id - patch_ring_osd_idx ") << std::endl;
    // for (int i = 0; i < patch_base_second_ring_verts.size(); ++i) {
    //     IndexType patch_ring_osd_idx = patch_base_second_ring_verts[i];
    //     IndexType patch_ring_kratos_id = r_patch_map_vert_to_node[0][patch_ring_osd_idx];
    //     KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: ") << i << " :          " << patch_ring_kratos_id << " - " << patch_ring_osd_idx << std::endl;
    // }
    // int n_patch_verts0 = patch_refiner->GetLevel(0).GetNumVertices();
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: num verts on level 0 of patch ") << n_patch_verts0 << std::endl;

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: loop over base_refined_second_ring_vertices ") << base_refined_second_ring_vertices << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: refined_base_ring_kratos_id - refined_base_ring_osd_idx ") << std::endl;
    // for (int i = 0; i < base_refined_second_ring_vertices.size(); ++i) {
    //     IndexType refined_base_ring_osd_idx = base_refined_second_ring_vertices[i];
    //     IndexType refined_base_ring_kratos_id = r_global_map_vert_to_node[1][refined_base_ring_osd_idx];
    //     KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: ") << i << " :          " << refined_base_ring_kratos_id << " - " << refined_base_ring_osd_idx << std::endl;
    // }
    // int n_vert_base_second_ring = base_refined_second_ring_vertices.size();
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: num verts in second ring of refined base ") << n_vert_base_second_ring << std::endl;
    // KRATOS_ERROR;

    // now patch_base_second_ring_verts contains all osd ids of the ctrl pts that influence the vertices of final regular patch
    // evaluate unit weight distribution for every node on the base patch level (refinement 0)
    std::map<IndexType, std::map<IndexType, double> > patch_refining_weights;
    GetRefiningWeights(patch_refining_weights, patch_base_second_ring_verts, regular_second_ring_vertices, r_patch_map_vert_to_node, patch_refiner, current_refining_depth);

    /// 2. get weights for original base control points of control polygon
    IndexType base_refined_osd_face_index = GetOSDIndexFromKratosID(FirstKratosFaceId, r_global_map_face_to_cond[1]);
    IndexType base_osd_face_index = refiner->GetLevel(1).GetFaceParentFace(base_refined_osd_face_index);
    IndexType kratos_control_polygon_base_face_index = r_global_map_face_to_cond[0][base_osd_face_index];
    // IndexType kratos_control_polygon_refined_base_face_index = r_global_map_face_to_cond[1][base_refined_osd_face_index];
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent base_osd_face_index on control polygon ") << base_osd_face_index << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: corresponding kratos_control_polygon_base_face_index ") << kratos_control_polygon_base_face_index << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FirstKratosFaceId  ") << FirstKratosFaceId << std::endl;
    std::vector<IndexType> base_second_ring_verts;    // contains osd ids
    std::vector<IndexType> base_refined_second_ring_verts;    // contains osd ids
    // base level
    GetSecondRingVertices(
        kratos_control_polygon_base_face_index, 
        base_second_ring_verts, 
        r_global_map_vert_to_node, 
        r_global_map_face_to_cond,
        refiner,
        0
    );
    // one refinement
    GetSecondRingVertices(
        FirstKratosFaceId, 
        base_refined_second_ring_verts, 
        r_global_map_vert_to_node, 
        r_global_map_face_to_cond,
        refiner,
        1
    );
    std::map<IndexType, std::map<IndexType, double> > base_refining_weights;
    GetRefiningWeights(base_refining_weights, base_second_ring_verts, base_refined_second_ring_verts, r_global_map_vert_to_node, refiner, 1);
    
    /// bring weights together
    const auto base_verts_begin = base_second_ring_verts.begin();
    const auto base_refined_verts_begin = base_refined_second_ring_verts.begin();
    const auto patch_verts_begin = patch_base_second_ring_verts.begin();
    const auto regular_patch_begin = regular_second_ring_vertices.begin();
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: combine weights ") << std::endl;
    for (IndexType i = 0; i < base_second_ring_verts.size(); ++i) {
        Vector weight_mult(16);
        weight_mult *= 0.0;
        auto base_v_osd_idx_it = base_verts_begin + i;
        IndexType base_kratos_node_id = r_global_map_vert_to_node[0][*base_v_osd_idx_it];
        // //
        // Vector w_vec = base_refining_weights[base_kratos_node_id];
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refining_weights.size ") << base_refining_weights.size() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_kratos_node_id ") << base_kratos_node_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refining_weights[base_kratos_node_id] ") << base_refining_weights[base_kratos_node_id] << std::endl;
        // // 
        for (IndexType k = 0; k < weight_mult.size(); ++k) {
            auto regular_patch_osd_idx_it = regular_patch_begin + k;
            IndexType regular_patch_kratos_node_id = r_patch_map_vert_to_node[current_refining_depth][*regular_patch_osd_idx_it];
            for (IndexType j = 0; j < patch_base_second_ring_verts.size(); ++j) {
                auto patch_v_osd_idx_it = patch_verts_begin + j;
                IndexType patch_base_kratos_node_id = r_patch_map_vert_to_node[0][patch_base_second_ring_verts[j]];
                IndexType refined_base_osd_node_id = GetOSDIndexFromKratosID(patch_base_kratos_node_id, r_global_map_vert_to_node[1]);
                IndexType refined_base_kratos_node_id = r_global_map_vert_to_node[1][base_refined_second_ring_verts[j]];
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: (patch_base_second_ring_verts[j] == *patch_v_osd_idx_it) ") << (patch_base_second_ring_verts[j] == *patch_v_osd_idx_it) << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_base_kratos_node_id ") << patch_base_kratos_node_id << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: refined_base_kratos_node_id ") << refined_base_kratos_node_id << std::endl;
                
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: indices -- i") << i << " -- j: " << j << " -- k: " << k << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refining_weights[base_kratos_node_id] ") << base_refining_weights[base_kratos_node_id] << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refining_weights[base_kratos_node_id][patch_base_kratos_node_id] ") << base_refining_weights[base_kratos_node_id][patch_base_kratos_node_id] << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_refining_weights.size ") << patch_refining_weights.size() << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_base_kratos_node_id ") << patch_base_kratos_node_id << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_refining_weights[patch_base_kratos_node_id] ") << patch_refining_weights[patch_base_kratos_node_id] << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_refining_weights[patch_base_kratos_node_id][regular_patch_kratos_node_id] ") << patch_refining_weights[patch_base_kratos_node_id][regular_patch_kratos_node_id] << std::endl;

                weight_mult[k] +=   base_refining_weights[base_kratos_node_id][patch_base_kratos_node_id] * 
                                    patch_refining_weights[patch_base_kratos_node_id][regular_patch_kratos_node_id];
                }
        }

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_kratos_node_id ") << base_kratos_node_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: weight_mult ") << weight_mult << std::endl;
        // weight_mult should have same ordering as final regular_second_ring_vertices
        RefiningWeights[base_kratos_node_id] = weight_mult;
    }    
    // KRATOS_INFO("RefiningWeights.size()") << RefiningWeights.size() << std::endl;
    
    // /// for debugging: compare weight_mult with result from refining original refiner (current_refining_depth+1)
    // RefineGeometry(rControlPolygon.GetMesh(), r_refined_patch_mp, r_global_map_vert_to_node, r_global_map_face_to_cond, current_refining_depth+1, refiner);
    // // skipping pushing to limit, not as exact as with pushing to limit
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: previous KratosFaceId ") << KratosFaceId << std::endl;
    // KratosFaceId = FindClosestConditionOnMesh(rFeNode, r_refined_patch_mp.GetMesh());
    // std::vector<IndexType> second_ring_verts;    // contains osd ids
    // bool is_irregular = GetSecondRingVertices(
    //     KratosFaceId, 
    //     second_ring_verts, 
    //     r_global_map_vert_to_node, 
    //     r_global_map_face_to_cond,
    //     refiner,
    //     current_refining_depth+1
    // );
    
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: debug KratosFaceId ") << KratosFaceId << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: first patch_is_irregular ") << patch_is_irregular << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: debug patch is_irregular ") << is_irregular << std::endl;
    // std::map< IndexType, std::map<IndexType, double> > debug_refining_weights;
    // GetRefiningWeights(debug_refining_weights, base_second_ring_verts, second_ring_verts, r_global_map_vert_to_node, refiner, current_refining_depth+1);
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_second_ring_verts ") << base_second_ring_verts << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_second_ring_verts ") << second_ring_verts << std::endl;
    // IndexType vector_idx = 0;
    // std::map<IndexType, Vector> DebugRefiningWeights;
    // for (auto const& vector_weights : debug_refining_weights) {
    //     // IndexType osd_base_ctrl_point_idx = base_second_ring_verts[vector_idx];
    //     // IndexType kratos_base_ctrl_point_id_1 = r_global_map_vert_to_node[0][osd_base_ctrl_point_idx];
    //     IndexType kratos_base_ctrl_point_id_2 = vector_weights.first;
    //     IndexType node_idx = 0;
    //     Vector weights(vector_weights.second.size());
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: debug_refining_weights loop number ") << vector_idx << std::endl;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: debug_refining_weights entry osd_base_ctrl_point_idx ") << osd_base_ctrl_point_idx << std::endl;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: debug_refining_weights entry osd_base_ctrl_point_idx ") << osd_base_ctrl_point_idx << std::endl;
    //     for (auto const& node_weight : vector_weights.second) {
    //         // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: node_idx ") << node_idx << std::endl;
    //         // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: node_weight.first ") << node_weight.first << std::endl;
    //         // Shouldn't the weigths be ordered like the numbering in second_ring_verts of the refined regular patch ?
    //         IndexType osd_patch_ctrl_point_idx = second_ring_verts[node_idx];
    //         IndexType kratos_patch_ctrl_point_id = r_global_map_vert_to_node[current_refining_depth+1][osd_patch_ctrl_point_idx];
    //         // weights[node_idx] = node_weight.second;
    //         weights[node_idx] = debug_refining_weights[kratos_base_ctrl_point_id_2][kratos_patch_ctrl_point_id];
    //         ++node_idx;
    //     }
    //     // RefiningWeights[kratos_base_ctrl_point_id_2] = weights;
    //     DebugRefiningWeights[kratos_base_ctrl_point_id_2] = weights;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: corresponding kratos_base_ctrl_point_id_1 ") << kratos_base_ctrl_point_id_1 << std::endl;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: corresponding kratos_base_ctrl_point_id ") << kratos_base_ctrl_point_id_2 << std::endl;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DebugRefiningWeights[kratos_base_ctrl_point_id] ") << weights << std::endl;
    //     // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: RefiningWeights[kratos_base_ctrl_point_id] ") << RefiningWeights[kratos_base_ctrl_point_id_2] << std::endl;
    //     ++vector_idx;
    // }
    // /// end debugging:

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rPatchNodesVector.size ") << rPatchNodesVector.size() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: end of function ") << std::endl;
    return rPatchNodesVector;
}


std::vector<NodeTypePointer> CreateRegularPatchForFeNode_Alternative(
    NodeType& rFeNode, 
    ModelPart& rControlPolygon,
    ModelPart& rGrid,
    ModelPart& rGridLimit,
    std::map<IndexType, Vector>& RefiningWeights,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_vert_to_node,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_face_to_cond,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode") << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rGrid.GetNode(1).Coordinates") << std::endl;
    // KRATOS_INFO(" -- ") << rGrid.GetNode(1).Coordinates()[0] << rGrid.GetNode(1).Coordinates()[1] << rGrid.GetNode(1).Coordinates()[2] << std::endl;
    SizeType max_refining_depth = 8;
    SizeType current_refining_depth = 0;    // we start at first refinement of control polygon
    Model patch_model;
    Model refined_patch_model;
    Model refined_limit_model;
    ModelPart& r_patch_mp = patch_model.CreateModelPart("patch", 3);
    ModelPart& r_refined_patch_mp = refined_patch_model.CreateModelPart("refined_patch", 3);
    ModelPart& r_refined_limit_mp = refined_limit_model.CreateModelPart("refined_limit", 3);
    r_patch_mp.CreateNewProperties(0);
    r_refined_patch_mp.CreateNewProperties(0);
    r_refined_limit_mp.CreateNewProperties(0);
    std::map<IndexType, std::vector<IndexType>> first_ring_vert_map;

    // 5.1 find closest control face on limit surface
    IndexType KratosFaceId = FindClosestConditionOnMesh(rFeNode, rGridLimit.GetMesh());
    IndexType FirstKratosFaceId = KratosFaceId;
    std::vector<IndexType> base_refined_second_ring_vertices;    // contains osd ids
    std::vector<IndexType> base_refined_second_ring_faces;       // contains osd ids
        
    // 5.2 create second ring vertices around control face
    bool patch_is_irregular = GetSecondRingVertices(
        KratosFaceId, 
        base_refined_second_ring_vertices, 
        r_global_map_vert_to_node,
        r_global_map_face_to_cond,
        refiner,
        1
    );
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refined_second_ring_vertices ") << base_refined_second_ring_vertices << std::endl;

    GetSecondRingFaces(
        KratosFaceId, 
        base_refined_second_ring_faces, 
        r_global_map_vert_to_node,
        r_global_map_face_to_cond,
        refiner,
        1
    );
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_refined_second_ring_faces ") << base_refined_second_ring_faces << std::endl;

    // set up local patch to be refined
    SetUpPatchMP(
        rGrid, 
        r_patch_mp, 
        base_refined_second_ring_vertices, 
        base_refined_second_ring_faces, 
        r_global_map_vert_to_node[1], 
        r_global_map_face_to_cond[1]
        );
    // r_patch_mp.SetNodes(rInputMesh.pNodes());
    // r_patch_mp.SetConditions(rInputMesh.pConditions());
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;
    // ModelPart& r_patch_mp = AuxiliarModelPartUtilities(rInputModelPart).DeepCopyModelPart("patch", &patch_model);
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: after DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;

    // NodeTypePointer as return type
    std::vector<NodeTypePointer> rPatchNodesVector;

    // in case first refinement was enough, fill rPatchNodesVector and RefiningWeights, then return
    // still need to adapt refinement weights for the case where only one refinement is needed
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: patch_is_irregular") << patch_is_irregular << std::endl;

    // optional: write files with refined patches
    std::string filename = r_patch_mp.Name() + "_ref" + std::to_string(current_refining_depth) + "_id" + std::to_string(rFeNode.GetId());
    r_patch_mp.Check();
    // OutputModelPartToVtk(r_patch_mp, filename);

    if (!patch_is_irregular) {

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: fill patch nodes vector with vertices from first refinement ") << std::endl;
        FillPatchNodesVector(rPatchNodesVector, rGrid.GetMesh(), base_refined_second_ring_vertices, r_global_map_vert_to_node[1]);

        // get weights from refining: meaning the influence of the base control points on the refined vertices of the regular patch
        IndexType osd_face_index = GetOSDIndexFromKratosID(KratosFaceId, r_global_map_face_to_cond[1]);
        OpenSubdiv::Far::TopologyLevel const & lastRefLevel = refiner->GetLevel(1);
        // get face index of control grid where FE node lies in
        for (SizeType ref = 0; ref < 1; ++ref) {
            osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
        }

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent osd_face_index on control polygon ") << osd_face_index << std::endl;
        IndexType kratos_base_face_index = r_global_map_face_to_cond[0][osd_face_index];
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent kratos_base_face_index on control polygon ") << kratos_base_face_index << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FirstKratosFaceId  ") << FirstKratosFaceId << std::endl;
        

        std::vector<IndexType> base_second_ring_vertices;    // contains osd ids
        // get second ring verts and faces for the base control grid P0
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: computing base_second_ring_vertices for refining weights due to 1 refinement ") << std::endl;
        GetSecondRingVertices(
            kratos_base_face_index, 
            base_second_ring_vertices, 
            r_global_map_vert_to_node, 
            r_global_map_face_to_cond,
            refiner,
            0
        );
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_second_ring_vertices ") << base_second_ring_vertices << std::endl;

        // now base_second_ring_vertices contains all osd ids of the ctrl pts that influence the vertices of final regular patch
        // evaluate unit weight distribution for every node on the base level (refinement 0)
        std::map<IndexType, std::map<IndexType, double> > refining_weights;
        GetRefiningWeights(refining_weights, base_second_ring_vertices, base_refined_second_ring_vertices, r_global_map_vert_to_node, refiner, 1);
        
        /// bring weights together
        const auto base_verts_begin = base_second_ring_vertices.begin();
        const auto base_refined_verts_begin = base_refined_second_ring_vertices.begin();
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: combine weights ") << std::endl;
        for (IndexType i = 0; i < base_second_ring_vertices.size(); ++i) {
            Vector weight_mult(16);
            weight_mult *= 0.0;
            auto base_v_osd_idx_it = base_verts_begin + i;
            IndexType base_kratos_node_id = r_global_map_vert_to_node[0][*base_v_osd_idx_it];

            for (IndexType k = 0; k < weight_mult.size(); ++k) {
                auto base_refined_verts_osd_idx_it = base_refined_verts_begin + k;
                IndexType regular_patch_kratos_node_id = r_global_map_vert_to_node[1][*base_refined_verts_osd_idx_it];
                weight_mult[k] += refining_weights[base_kratos_node_id][regular_patch_kratos_node_id];
            }

            RefiningWeights[base_kratos_node_id] = weight_mult;
        }

        return rPatchNodesVector;
    }
    // else if we have irregular patch:


    // refine the patch until fe-node is inside regular patch
    std::map<IndexType, std::map<IndexType, IndexType>> r_patch_map_vert_to_node;
    std::map<IndexType, std::map<IndexType, IndexType>> r_patch_map_face_to_cond;
    OpenSubdiv::Far::TopologyRefiner* patch_refiner = nullptr;

    std::vector<IndexType> regular_second_ring_vertices;    // regular patch
    while (patch_is_irregular && max_refining_depth > current_refining_depth) 
    {
        patch_refiner = GenerateOpenSubdivRefiner(r_patch_mp.GetMesh(), r_patch_map_vert_to_node, r_patch_map_face_to_cond, false);
        current_refining_depth += 1;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: current_refining_depth of patch") << current_refining_depth << std::endl;

        // auto geom = rInputMesh.GetCondition(KratosFaceId).GetGeometry();
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom ") << geom << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom IDs ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;
        RefineGeometry(r_patch_mp.GetMesh(), r_refined_patch_mp, r_patch_map_vert_to_node, r_patch_map_face_to_cond, 1, patch_refiner);
        // TODO: compute new weights automatically when refining

        std::map<IndexType, std::vector<IndexType>> first_ring_vert_map;
        GetFirstRingVertices(r_refined_patch_mp.GetMesh(), first_ring_vert_map, r_patch_map_vert_to_node, r_patch_map_face_to_cond, patch_refiner, 1);
        // create modelpart whose nodes lie on the limit surface
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << r_refined_patch_mp.Name() << std::endl;
        ModelPart& r_refined_limit_mp = AuxiliarModelPartUtilities(r_refined_patch_mp).DeepCopyModelPart(r_refined_patch_mp.Name()+"_limit_"+std::to_string(current_refining_depth), &refined_limit_model);
        PushNodesToLimit(r_refined_patch_mp, r_refined_limit_mp, first_ring_vert_map);
        
        // after refining: check which face fe_node lies in
        KratosFaceId = FindClosestConditionOnMesh(rFeNode, r_refined_limit_mp.GetMesh());
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New KratosFaceId in refined mesh ") << KratosFaceId << std::endl;

        // create 2-ring vertex neighbourhoods around the face the fe-node lies in
        patch_is_irregular = GetSecondRingVertices(
            KratosFaceId, 
            regular_second_ring_vertices, 
            r_patch_map_vert_to_node,
            r_patch_map_face_to_cond,
            patch_refiner,
            1
            );
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New patch for KratosFaceId: patch_is_irregular ") << patch_is_irregular << std::endl;
        // optional: write files with refined patches
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: r_refined_patch_mp.NumberOfNodes() ") << r_refined_patch_mp.NumberOfNodes() << std::endl;
        std::string filename = r_refined_patch_mp.Name() + "_ref" + std::to_string(current_refining_depth) + "_id" + std::to_string(rFeNode.GetId());
        r_refined_patch_mp.Check();
        // OutputModelPartToVtk(r_refined_patch_mp, filename);

        // second ring faces to create new patch model part
        std::vector<IndexType> refined_second_ring_faces;       // contains osd ids
        GetSecondRingFaces(
            KratosFaceId, 
            refined_second_ring_faces, 
            r_patch_map_vert_to_node,
            r_patch_map_face_to_cond,
            refiner,
            1
        );
        // set up local patch to be refined
        r_patch_mp.Clear();
        SetUpPatchMP(
            rGrid, 
            r_patch_mp, 
            regular_second_ring_vertices, 
            refined_second_ring_faces, 
            r_patch_map_vert_to_node[1], 
            r_patch_map_face_to_cond[1]
            );
    }

    if (current_refining_depth >= max_refining_depth && patch_is_irregular) {
        // throw error
        KRATOS_ERROR << "Reached maximum refining depth, could not determine regular patch for FE-node " << rFeNode.GetId() << std::endl;

        // actually near enoug to an irregular vertex to just take limit map of the neighbourhood
        // 1. get nearest node on irregular patch to fe-node (should be irregular vertex)
        // 2. evaluate limit map in 1-ring neighb of irregular vertex to get weights of surrounding vertices
    }
    // else regular patch was found after refining 
    

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: fill patch nodes vector with vertices from refinement level ") << current_refining_depth << std::endl;
    FillPatchNodesVector(rPatchNodesVector, r_refined_patch_mp.GetMesh(), regular_second_ring_vertices, r_patch_map_vert_to_node[current_refining_depth]);

    // get weights from refining: meaning the influence of the base control points on the refined vertices of the regular patch
    /// 1. get weights for base level of patch_refiner
    IndexType osd_face_index = GetOSDIndexFromKratosID(KratosFaceId, r_patch_map_face_to_cond[current_refining_depth]);
    OpenSubdiv::Far::TopologyLevel const & lastRefLevel = patch_refiner->GetLevel(current_refining_depth);
    // get face index of control grid where FE node lies in
    for (SizeType ref = 0; ref < current_refining_depth; ++ref) {
        osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
    }

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent osd_face_index on control polygon ") << osd_face_index << std::endl;
    IndexType kratos_patch_base_face_index = r_patch_map_face_to_cond[0][osd_face_index];

    std::vector<IndexType> patch_base_second_ring_verts;    // contains osd ids
    // get second ring verts and faces for the base control grid P0
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: compute patch_base_second_ring_verts to get refining weights for refinement number ") << current_refining_depth << std::endl;
    GetSecondRingVertices(
        FirstKratosFaceId, 
        patch_base_second_ring_verts, 
        r_patch_map_vert_to_node, 
        r_patch_map_face_to_cond,
        patch_refiner,
        0
    );
    
    // now patch_base_second_ring_verts contains all osd ids of the ctrl pts that influence the vertices of final regular patch
    // evaluate unit weight distribution for every node on the base patch level (refinement 0)
    std::map<IndexType, std::map<IndexType, double> > patch_refining_weights;
    GetRefiningWeights(patch_refining_weights, patch_base_second_ring_verts, regular_second_ring_vertices, r_patch_map_vert_to_node, patch_refiner, current_refining_depth);

    /// 2. get weights for original base control points of control polygon
    IndexType base_refined_osd_face_index = GetOSDIndexFromKratosID(FirstKratosFaceId, r_global_map_face_to_cond[1]);
    IndexType base_osd_face_index = refiner->GetLevel(1).GetFaceParentFace(base_refined_osd_face_index);
    IndexType kratos_control_polygon_base_face_index = r_global_map_face_to_cond[0][base_osd_face_index];
    // IndexType kratos_control_polygon_refined_base_face_index = r_global_map_face_to_cond[1][base_refined_osd_face_index];
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent base_osd_face_index on control polygon ") << base_osd_face_index << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: corresponding kratos_control_polygon_base_face_index ") << kratos_control_polygon_base_face_index << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: FirstKratosFaceId  ") << FirstKratosFaceId << std::endl;
    std::vector<IndexType> base_second_ring_verts;    // contains osd ids
    std::vector<IndexType> base_refined_second_ring_verts;    // contains osd ids
    // base level
    GetSecondRingVertices(
        kratos_control_polygon_base_face_index, 
        base_second_ring_verts, 
        r_global_map_vert_to_node, 
        r_global_map_face_to_cond,
        refiner,
        0
    );
    // one refinement
    GetSecondRingVertices(
        FirstKratosFaceId, 
        base_refined_second_ring_verts, 
        r_global_map_vert_to_node, 
        r_global_map_face_to_cond,
        refiner,
        1
    );
    std::map<IndexType, std::map<IndexType, double> > base_refining_weights;
    GetRefiningWeights(base_refining_weights, base_second_ring_verts, base_refined_second_ring_verts, r_global_map_vert_to_node, refiner, 1);
    
    /// bring weights together
    const auto base_verts_begin = base_second_ring_verts.begin();
    const auto base_refined_verts_begin = base_refined_second_ring_verts.begin();
    const auto patch_verts_begin = patch_base_second_ring_verts.begin();
    const auto regular_patch_begin = regular_second_ring_vertices.begin();
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: combine weights ") << std::endl;
    for (IndexType i = 0; i < base_second_ring_verts.size(); ++i) {
        Vector weight_mult(16);
        weight_mult *= 0.0;
        auto base_v_osd_idx_it = base_verts_begin + i;
        IndexType base_kratos_node_id = r_global_map_vert_to_node[0][*base_v_osd_idx_it];

        for (IndexType k = 0; k < weight_mult.size(); ++k) {
            auto regular_patch_osd_idx_it = regular_patch_begin + k;
            IndexType regular_patch_kratos_node_id = r_patch_map_vert_to_node[current_refining_depth][*regular_patch_osd_idx_it];
            for (IndexType j = 0; j < patch_base_second_ring_verts.size(); ++j) {
                auto patch_v_osd_idx_it = patch_verts_begin + j;
                IndexType patch_base_kratos_node_id = r_patch_map_vert_to_node[0][patch_base_second_ring_verts[j]];
                IndexType refined_base_osd_node_id = GetOSDIndexFromKratosID(patch_base_kratos_node_id, r_global_map_vert_to_node[1]);
                IndexType refined_base_kratos_node_id = r_global_map_vert_to_node[1][base_refined_second_ring_verts[j]];

                weight_mult[k] +=   base_refining_weights[base_kratos_node_id][patch_base_kratos_node_id] * 
                                    patch_refining_weights[patch_base_kratos_node_id][regular_patch_kratos_node_id];
                }
        }

        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_kratos_node_id ") << base_kratos_node_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: weight_mult ") << weight_mult << std::endl;
        // weight_mult should have same ordering as final regular_second_ring_vertices
        RefiningWeights[base_kratos_node_id] = weight_mult;


    }    

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rPatchNodesVector.size ") << rPatchNodesVector.size() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: end of function ") << std::endl;
    return rPatchNodesVector;
}



CoordinatesArrayType ComputeGlobalFeCoordsBasedOnWeights(Vector& rWeights, PointerVector<NodeType> pPoints)
{
    CoordinatesArrayType global_coords;
    global_coords[0] = 0.0;
    global_coords[1] = 0.0;
    global_coords[2] = 0.0;
    for (SizeType i = 0; i < rWeights.size(); ++i) {
        global_coords[0] += rWeights[i] * pPoints[i].Coordinates()[0];
        global_coords[1] += rWeights[i] * pPoints[i].Coordinates()[1];
        global_coords[2] += rWeights[i] * pPoints[i].Coordinates()[2];
    }
    return global_coords;
}

void CatmullClarkSDS::CreateMappingMatrix(
    std::vector<double>& rOutputData,
    ModelPart& rControlPolygon,
    const ModelPart& rControlledMesh,
    const bool FixFreeEdges) 
{
    // set output vector containing mapping weights to correct size
    rOutputData.resize(rControlPolygon.NumberOfNodes() * rControlledMesh.NumberOfNodes());
    // create boundary ModelPart for open geometries
    // for now we omit creating boundary edges
    // ModelPart& rBoundaryModelPart = ModelPart(rControlPolygon.Name() + "_edges", &rControlPolygon.mpVariablesList, rControlPolygon);       // rControlPolygon.HasSubModelPart(rControlPolygon.Name() + "_edges") ? rControlPolygon.GetSubModelPart(rControlPolygon.Name() + "_edges") : rControlPolygon.CreateSubModelPart(rControlPolygon.Name() + "_edges");
    // ExtractEdgeNodes(rControlPolygon.Name() + "_edges");

    // Do one subdivision step to work only on quads, bcs for now the implementation is hardcoded for quads
    // ModelPart::MeshType r_first_subd_for_quads;

    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: original number of Control Nodes ") << rControlPolygon.NumberOfNodes() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: original number of Control Conds ") << rControlPolygon.NumberOfConditions() << std::endl;


    Model model;
    Model model_limit;
    ModelPart& r_first_subd_for_quads = model.CreateModelPart("first_subd", 3);
    r_first_subd_for_quads.CreateNewProperties(0);

    // the maps store for every refinement level the relation between Kratos IDs and OSD indices
    // structure is the following:
    // map<level, map<KratosID, OSDIndex>>
    std::map<IndexType, std::map<IndexType, IndexType>> global_map_osd_vert_to_kratos_node_ids;
    std::map<IndexType, std::map<IndexType, IndexType>> global_map_osd_face_to_kratos_cond_ids;

    /// 1. create main opensubdiv refiner for generation of first and 
    /// second ring vertices of base control point grid
    OpenSubdiv::Far::TopologyRefiner * refiner_main = GenerateOpenSubdivRefiner(rControlPolygon.GetMesh(), global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, FixFreeEdges);

    // 2. refine control polygon once in order to have only quads
    // main reasons are that pushing nodes to limit currently only works with quads
    RefineGeometry(rControlPolygon.GetMesh(), r_first_subd_for_quads, global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, 1, refiner_main);

    /// debugging
    std::string filename = r_first_subd_for_quads.Name();
    r_first_subd_for_quads.Check();
    // OutputModelPartToVtk(r_first_subd_for_quads, filename);
    /// debugging

    // 3. create first ring vertices around every control node for pushing to limit 
    // (done once for entore control polygon before loop over fe nodes)
    std::map<IndexType, std::vector<IndexType>> first_ring_vertices;
    // std::map<IndexType, std::vector<IndexType>> base_cp_first_ring_neighb, base_cp_second_ring_neighb;
    GetFirstRingVertices(r_first_subd_for_quads.GetMesh(), first_ring_vertices, global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, refiner_main, 1);

    // 4. push control nodes to limit to enable search for control faces that fe-nodes lie in
    ModelPart& r_points_on_limit = AuxiliarModelPartUtilities(r_first_subd_for_quads).DeepCopyModelPart(r_first_subd_for_quads.Name()+"_limit", &model_limit);
    PushNodesToLimit(r_first_subd_for_quads, r_points_on_limit, first_ring_vertices);

    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rControlPolygon num nodes ") << rControlPolygon.NumberOfNodes() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: total number of FE Nodes ") << rControlledMesh.NumberOfNodes() << std::endl;

    // 5. start loop over FE-nodes in FE-mesh
    const auto fe_node_begin = rControlledMesh.NodesBegin();
// #pragma omp parallel for
    for (IndexType fe_node_index = 0; fe_node_index < rControlledMesh.NumberOfNodes(); ++fe_node_index) {
        IndexType row_index = fe_node_index * rControlPolygon.NumberOfNodes();
        auto fe_node_it = fe_node_begin + fe_node_index;
        // NodeType& fe_node = *fe_node_it;
        IndexType fe_node_id = fe_node_it->GetId();
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: starting with fe node id ") << fe_node_id << std::endl;

        // full section moved inside CreateRegularPatchForFeNode
        // // 5.1 find closest control face on limit surface
        // IndexType kratos_face_id = FindClosestConditionOnMesh(*fe_node_it, r_points_on_limit.GetMesh());

        // // 5.2 create second ring vertices around control face
        // IndexType current_ref = 1;
        // std::vector<IndexType> second_ring_vertices;    // contains osd ids
        // bool irregular_patch = GetSecondRingVertices(
        //     kratos_face_id,
        //     second_ring_vertices,
        //     global_map_osd_vert_to_kratos_node_ids,
        //     global_map_osd_face_to_kratos_cond_ids,
        //     refiner_main,
        //     current_ref
        // );

        std::vector<IndexType> patch_ids;
        std::vector<NodeTypePointer> patch;
        
                    // create refining weight values for refined regular patches if patch is irregular
                    std::map<IndexType, Vector> refining_weights_map;

                    // create regular patch by refinement
                    patch = CreateRegularPatchForFeNode(
                        // fe_node, 
                        *fe_node_it, 
                        rControlPolygon,
                        r_first_subd_for_quads, 
                        r_points_on_limit,
                        refining_weights_map, 
                        global_map_osd_vert_to_kratos_node_ids, 
                        global_map_osd_face_to_kratos_cond_ids,
                        refiner_main
                    );
                        // r_first_subd_for_quads, 
                    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: returned patch ") << patch << std::endl;

                    // nur vorübergehend, um r_first_subd_for_quads zu füllen
                    // RefineGeometry(rControlPolygon.GetMesh(), r_first_subd_for_quads, FixFreeEdges);

                    // special case where we found the fe node to be extremely near to an irregular point
                    // KRATOS_INFO("patch size is") << patch.size() << std::endl;
                    if (patch.size() != 16) {
                        // KRATOS_INFO("first entry of refining_weights_map size is") << refining_weights_map[0].size() << std::endl;
                        // frst check if all dimensions coincide
                        KRATOS_ERROR_IF_NOT(patch.size() == refining_weights_map[0].size());
                        // else fill mapping matrix with weights for corresponding control points
                        Vector weight_values = refining_weights_map[0];
                        // KRATOS_INFO("first weight_values") << weight_values << std::endl;
                        
                        const auto& patch_begin = patch.begin();
                        for (IndexType i = 0; i < patch.size(); ++i) {
                            const auto& patch_node_it = patch_begin + i;
                            IndexType ctrl_node_id = (*patch_node_it)->GetId();
                            rOutputData[row_index + ctrl_node_id-1] = weight_values[i];
                        }
                        // KRATOS_INFO("rOutputData") << rOutputData << std::endl;

                        // KRATOS_INFO("Skipping building regular patch, since fe-node is very near to irregular vertex!") << std::endl;
                        // got to next fe node
                        // KRATOS_ERROR << "Error out to check out data" << std::endl;

                        continue;
                    }
                    // else patch.size() == 16 and we build the regular patch 
                    // to find fe node's parametric coordinates
        

        // 2. then find the parametric coordinates of the euclidean coordinates on limit surface with Newton-Raphson
        //      2.1 build b-spline surface with regular patch
                // Create bivariate nurbs surface.
                typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
                PointerVector<NodeType> points; // fill points with regular patch
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: patch.size() ") << patch.size() << std::endl;
                for(IndexType j = 0; j < patch.size(); j++) {
                    NodeTypePointer p_node_j = patch[j];
                    points.push_back(Kratos::make_intrusive<NodeType>(j, p_node_j->X(), p_node_j->Y(), p_node_j->Z() ));
                    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: p_node_j ") << p_node_j << std::endl;
                }

                SizeType OrderU = 3;
                SizeType OrderV = 3;
                SizeType num_cp_u = 4;
                SizeType num_cp_v = 4;
                SizeType number_of_knots_u = OrderU + num_cp_u + 1;
                SizeType number_of_knots_v = OrderV + num_cp_v + 1;
                Vector knot_vector_u(number_of_knots_u);
                Vector knot_vector_v(number_of_knots_v);

                double incr = 1.0/(number_of_knots_u-1.0);

                for(SizeType j = 0; j < number_of_knots_u; j++) {
                    knot_vector_u[j] = double(j)*incr;
                    knot_vector_v[j] = double(j)*incr;
                }
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: knot_vector_u ") << knot_vector_u << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: knot_vector_v ") << knot_vector_v << std::endl;


                auto p_surface_geometry = Kratos::make_shared<NurbsSurfaceGeometryType>(
                    points, OrderU, OrderV, knot_vector_u, knot_vector_v);

        //      2.2 solve for parametric coordinates
                CoordinatesArrayType local_coords;
                p_surface_geometry->ProjectionPointGlobalToLocalSpace(fe_node_it->Coordinates(), local_coords);
                // if (limit_intersection_points_map.find(fe_node.GetId()) == limit_intersection_points_map.end()) {
                //     p_surface_geometry->ProjectionPointGlobalToLocalSpace(fe_node.Coordinates(), local_coords);
                // }
                // else {
                //     CoordinatesArrayType intersection_coords = limit_intersection_points_map[fe_node.GetId()].Coordinates();
                //     p_surface_geometry->ProjectionPointGlobalToLocalSpace(intersection_coords, local_coords);
                // }
                // useful debug:
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: FE global_coords   ") << fe_node_it->Coordinates() << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: FE local_coords    ") << local_coords << std::endl;

        // 3. evaluate the fe point on limit surface to establish relation to (weights of) the control points
                Vector patch_weights;
                p_surface_geometry->ShapeFunctionsValues(patch_weights, local_coords);
                
                Vector refining_weights(refining_weights_map.size());
                PointerVector<NodeType> base_points;
                IndexType index = 0;
                for (auto const& ref_w_it : refining_weights_map) {

                    refining_weights[index] = inner_prod(patch_weights, ref_w_it.second);

                    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: ref_w_it.first   ") << ref_w_it.first << std::endl;
                    NodeTypePointer p_node = rControlPolygon.pGetNode(ref_w_it.first);
                    base_points.push_back( Kratos::make_intrusive<NodeType>(index, p_node->X(), p_node->Y(), p_node->Z() ));

                    ++index;
                }

                CoordinatesArrayType global_coords_nurbs = ComputeGlobalFeCoordsBasedOnWeights(patch_weights, points);
                // useful debug:
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: recomputed global_coords_nurbs with shape functions") << global_coords_nurbs << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: patch_weights ") << patch_weights << std::endl;
                CoordinatesArrayType global_coords_ref_weights = ComputeGlobalFeCoordsBasedOnWeights(refining_weights, base_points);
                // useful debug:
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: recomputed global_coords_ref_weights with shape functions") << global_coords_ref_weights << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: refining_weights ") << refining_weights << std::endl;
                // CoordinatesArrayType difference;
                // difference = fe_node_it->Coordinates() - global_coords_nurbs;
                // double diff_norm = difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2];
                // diff_norm = sqrt(diff_norm);
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: diff_norm ") << diff_norm << std::endl;
                // if (diff_norm > 0.5) {
                //     KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: control points ") << points << std::endl;
                // }

        //      3.1 save the relation in mapping matrix
                // todo : treat refined regular patches

                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: refining_weights_map.size() ") << refining_weights_map.size() << std::endl;

                for (auto const& ctrl_pt_it : refining_weights_map) {
                    // KRATOS_INFO(" ") << ctrl_pt_it.first << " : " << ctrl_pt_it.second << std::endl;

                    IndexType ctrl_node_id = ctrl_pt_it.first;
                    // KRATOS_INFO(" row_index ") << row_index << std::endl;
                    // KRATOS_INFO(" control_node_id ") << ctrl_node_id << std::endl;
                    // KRATOS_INFO(" entry index = row_index + control_node_id-1 ") << row_index + ctrl_node_id-1 << std::endl;
                    // KRATOS_INFO(" inner_prod(patch_weights, ctrl_pt_it.second) ") << inner_prod(patch_weights, ctrl_pt_it.second) << std::endl;
                    // output should have data structure: matrix dim = [num_fe_nodes x num_ctrl_nodes]
                    // thus every row takes num_ctrl_pts entries
                    // since we do outer loop over fe nodes, we need to set entries for the column here
                    rOutputData[row_index + ctrl_node_id-1] = inner_prod(patch_weights, ctrl_pt_it.second);
                }


                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rOutputData ") << rOutputData << std::endl;
                
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: fe_node_index ") << fe_node_index << std::endl;

                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: row_index ") << row_index << std::endl;


                bool error_out_in_first_loop = false;
                KRATOS_ERROR_IF(error_out_in_first_loop);
                // KRATOS_ERROR_IF(fe_node_index == 0);

        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: successfully finished loop number ") << fe_node_it->GetId() << std::endl;

    }
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rOutputData ") << rOutputData << std::endl;
    KRATOS_INFO("Finished creating mapping matrix...") << std::endl;
};

} // namespace kratos