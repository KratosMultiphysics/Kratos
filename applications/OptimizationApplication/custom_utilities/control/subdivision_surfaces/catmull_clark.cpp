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
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
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

using IndexType = std::size_t;
using NodeType = Node::NodeType;
using GeometryType = Geometry<Node>;
using ModelPart = ModelPart;
// using MeshType = ModelPart::MeshType;
using CoordinatesArrayType = Point::CoordinatesArrayType;
typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

// typedefs for nearest point search
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

void CreateFirstRingNeighbourhoodsOfNodes(
    ModelPart::MeshType& rInputMesh,
    std::map<IndexType, std::vector<IndexType>>& rOutputMap,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_verts,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_faces,
    OpenSubdiv::Far::TopologyRefiner* refiner = nullptr
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
    KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes") << std::endl;
    std::map<IndexType, std::vector<IndexType>> first_ring_node_neighbourhoods;

    SizeType ref_level = 0;
    // if refiner not passed, generate a new one based on rInputMesh
    if (!refiner) {
        // std::map<IndexType, IndexType> map_verts, map_conds;
        refiner = GenerateOpenSubdivRefiner(rInputMesh, r_global_map_verts, r_global_map_faces, false);
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: no refiner specified... creating new one") << std::endl;
    }
    else { 
        ref_level = refiner->GetNumLevels() - 1; 
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: using specified refiner") << std::endl;
    }

    OpenSubdiv::Far::TopologyLevel const & refLevel = refiner->GetLevel(ref_level);
    // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: current refiner level ") << ref_level << std::endl;

    // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: rInputMesh.NumberOfNodes() ") << rInputMesh.NumberOfNodes() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: refLevel.GetNumVertices() ") << refLevel.GetNumVertices() << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: refLevel.GetNumFaces() ") << refLevel.GetNumFaces() << std::endl;
    SizeType num_vertices = refLevel.GetNumVertices();
    SizeType num_nodes = rInputMesh.NumberOfNodes();
    KRATOS_ERROR_IF_NOT(num_nodes == num_vertices);

    std::map<IndexType, IndexType> map_verts = r_global_map_verts[ref_level];
    std::map<IndexType, IndexType> map_faces = r_global_map_faces[ref_level];

    // initializes the 1-ring node neighbourhoods for every node in the input mesh
    
    // loop over all nodes to find neighbourhood of every node
    IndexType vert_index = 0;
    for (auto node_it = rInputMesh.NodesBegin(); node_it != rInputMesh.NodesEnd(); ++node_it) {
        // find control node in mesh rInputMesh
        const NodeType& r_center_node = *node_it;
        IndexType kratos_center_node_id = r_center_node.GetId();
        IndexType osd_center_vert_idx = GetOSDIndexFromKratosID(kratos_center_node_id, map_verts);
        // IndexType osd_center_vertex_index = r_center_node_id - 1;
        // skip free edge nodes since they need to be evaluated seperately
        if (refLevel.IsVertexBoundary(osd_center_vert_idx)) {
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
                // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: debug 5 : *missing_it ") << *missing_it << std::endl;
            }
        }
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: temp_neighbourhood.size() ") << temp_neighbourhood.size() << std::endl;
        rOutputMap[kratos_center_node_id] = temp_neighbourhood;
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: rOutputMap[r_center_node_id] ") << rOutputMap[r_center_node_id] << std::endl;
    }
    // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfNodes :: rOutputMap.size() ") << rOutputMap.size() << std::endl;
    // for(const auto& elem : rOutputMap)
    // {
    //     KRATOS_INFO("CATMULL_CLARK                 :: ") << elem.first << " : " << elem.second << std::endl;
    // }
};

void CreateModelPartFor2RingPatch(
    ModelPart::MeshType& rInputMesh, 
    ModelPart& rPatchModelPart, 
    std::vector<IndexType> SecondRingVerts, // contains osd ids
    std::vector<IndexType> SecondRingFaces  // contains osd ids
    ) 
{
    // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch") << std::endl;

    // create nodes
    NodesContainerType::Pointer p_new_nodes_container = make_shared<NodesContainerType>();
    for (auto vert_it = SecondRingVerts.begin(); vert_it != SecondRingVerts.end(); ++vert_it) {
        IndexType vertex_osd_index = *vert_it;
        IndexType kratos_node_index = vertex_osd_index + 1;
        std::vector<double> node_coord(3);
        node_coord[0] = rInputMesh.GetNode(kratos_node_index).Coordinates()[0];
        node_coord[1] = rInputMesh.GetNode(kratos_node_index).Coordinates()[1];
        node_coord[2] = rInputMesh.GetNode(kratos_node_index).Coordinates()[2];
        NodeTypePointer p_node = make_intrusive<Node>(kratos_node_index, node_coord);
        // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: kratos_node_index") << kratos_node_index << std::endl;
        p_new_nodes_container->push_back(p_node);
    }
    rPatchModelPart.SetNodes(p_new_nodes_container);
    // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: rPatchModelPart.NumberOfNodes()") << rPatchModelPart.NumberOfNodes() << std::endl;

    // create conditions
    std::string condition_name = "SurfaceCondition3D4N";
    for (auto face_it = SecondRingFaces.begin(); face_it != SecondRingFaces.end(); ++face_it) {
        IndexType face_osd_index = *face_it;
        IndexType kratos_cond_index = face_osd_index + 1;
        // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: kratos_cond_index") << kratos_cond_index << std::endl;
        ConditionType& cond = rInputMesh.GetCondition(kratos_cond_index);
        SizeType number_of_nodes = cond.GetGeometry().size();
        // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: number_of_nodes") << number_of_nodes << std::endl;
        std::vector<IndexType> node_indices(number_of_nodes);
        for (int i = 0; i < number_of_nodes; ++i) {
            node_indices[i] = cond.GetGeometry()[i].GetId();
        }
        // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: node_indices") << node_indices << std::endl;

        Properties::Pointer p_property = rPatchModelPart.pGetProperties(0);
        rPatchModelPart.CreateNewCondition(condition_name, kratos_cond_index, node_indices, p_property);
    }
    // KRATOS_INFO("CATMULL_CLARK :: CreateModelPartFor2RingPatch :: rPatchModelPart.NumberOfConditions()") << rPatchModelPart.NumberOfConditions() << std::endl;
}

bool GetSecondRingVerticesAndFaces(
    IndexType KratosFaceId, 
    ModelPart::MeshType& rInputMesh,
    std::vector<IndexType>& SecondRingVertices, // contains osd ids
    std::vector<IndexType>& SecondRingFaces     // contains osd ids
    )
{
    SecondRingVertices.clear();
    SecondRingFaces.clear();
    KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces") << std::endl;
    // assuming only quadrilaterals, since at least 1 subdivision has been performed previously
    bool is_irregular_patch = false;

    std::map<IndexType, std::map<IndexType, IndexType>> map_vert_to_node, map_face_to_cond;
    OpenSubdiv::Far::TopologyRefiner * refiner_for_geometry_description = GenerateOpenSubdivRefiner(rInputMesh, map_vert_to_node, map_face_to_cond, false);
    
    IndexType osd_face_id = GetOSDIndexFromKratosID(KratosFaceId, map_face_to_cond[0]);
    KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: KratosFaceId") << KratosFaceId << std::endl;
    KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: osd_face_id") << osd_face_id << std::endl;

    // get local patch geometry around osd_face_id
    // std::vector<int> num_vertices_per_face;
    // std::vector<int> vertex_indices_per_face;

    OpenSubdiv::Far::TopologyLevel const & firstLevel = refiner_for_geometry_description->GetLevel(0);

    // second ring vertices and faces
    SecondRingFaces.push_back(osd_face_id);
    ConstIndexArray vert_indices_of_face = firstLevel.GetFaceVertices(osd_face_id);
    // vert_indices_of_face.size() should be 4 since there has been at least 1 subdivision previously

    // check if irregular
    for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
        if (!firstLevel.IsVertexValenceRegular(vert_indices_of_face[vert_i])) is_irregular_patch = true;
        
        ConstIndexArray faces_around_vert = firstLevel.GetVertexFaces(vert_indices_of_face[vert_i]);
        for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                ConstIndexArray vertices_of_face_j = firstLevel.GetFaceVertices(face_index);
                if (vertices_of_face_j.size() != 4) is_irregular_patch = true;
        }
    }

    // create 2-ring neighbourhoods
    if (is_irregular_patch) {
        KRATOS_INFO("GetSecondRingVerticesAndFaces :: irregular patch ...") << std::endl;

        for (IndexType vert_i = 0; vert_i < vert_indices_of_face.size(); ++vert_i) {
            IndexType vert_index = vert_indices_of_face[vert_i];
            ConstIndexArray faces_around_vert = firstLevel.GetVertexFaces(vert_index);

            for (IndexType face_j = 0; face_j < faces_around_vert.size(); ++face_j) {
                IndexType face_index = faces_around_vert[face_j];
                // skip if face already accounted for, else add to vector
                auto f_it = std::find(SecondRingFaces.begin(), SecondRingFaces.end(), face_index);
                if (f_it == SecondRingFaces.end()) SecondRingFaces.push_back(face_index);
                // num_vertices_per_face.push_back(vertices_of_face_j.size());

                ConstIndexArray vertices_of_face_j = firstLevel.GetFaceVertices(face_index);
                // KRATOS_INFO("Opensubdiv vertices of face ") << face_index << "\n [ " << vertices_of_face_j[0] << " , " << vertices_of_face_j[1] << " , " << vertices_of_face_j[2] << " , " << vertices_of_face_j[3] << " ]" << std::endl;
                // auto cond_geometry = rInputMesh.GetCondition(face_index+1).GetGeometry();
                // KRATOS_INFO("Kratos Geometry for condition ") << face_index+1 << "\n [ " << cond_geometry[0].GetId() << " , " << cond_geometry[1].GetId() << " , " << cond_geometry[2].GetId() << " , " << cond_geometry[3].GetId() << " ]" << std::endl;
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
        KRATOS_INFO("GetSecondRingVerticesAndFaces :: regular patch ...") << std::endl;
        // prepare vector structure for regular patch
        SecondRingVertices.resize(16);

        std::vector<std::vector<IndexType>> patchPointsIndicesPerFace = { {  5,  4,  0,  1 },
                                                                        {  6,  2,  3,  7 },
                                                                        { 10, 11, 15, 14 },
                                                                        {  9, 13, 12,  8 } };
        
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: Found regular patch, creating ordered patch... ") << std::endl;

        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: current osd_face_id ") << osd_face_id << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: vertices of osd_face_id ") << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [0] ") << firstLevel.GetFaceVertices(osd_face_id)[0] << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [1] ") << firstLevel.GetFaceVertices(osd_face_id)[1] << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [2] ") << firstLevel.GetFaceVertices(osd_face_id)[2] << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [3] ") << firstLevel.GetFaceVertices(osd_face_id)[3] << std::endl;

        // create ordered regular patch
        IndexType patch_position_index = 0;
        for (auto vert_index_it = vert_indices_of_face.begin(); vert_index_it != vert_indices_of_face.end(); ++vert_index_it) {
            IndexType vert_index = *vert_index_it;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: Current vert_index ") << vert_index << std::endl;
            std::vector<IndexType> patch_indices = patchPointsIndicesPerFace[patch_position_index];
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: patch_position_index") << patch_position_index << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: patch_indices") << patch_indices << std::endl;
            patch_position_index += 1;

            ConstIndexArray face_indices = firstLevel.GetVertexFaces(vert_index);
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: face_indices for vert_index") << vert_index << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [0] ") << face_indices[0] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [1] ") << face_indices[1] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [2] ") << face_indices[2] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: [3] ") << face_indices[3] << std::endl;
            IndexType opposite_face_index = GetOppositeIndex(face_indices, osd_face_id);
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: opposite_face_index") << opposite_face_index << std::endl;
            ConstIndexArray vert_indices_for_patch = firstLevel.GetFaceVertices(opposite_face_index);
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: vert_indices_for_patch[0]") << vert_indices_for_patch[0] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: vert_indices_for_patch[1]") << vert_indices_for_patch[1] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: vert_indices_for_patch[2]") << vert_indices_for_patch[2] << std::endl;
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: vert_indices_for_patch[3]") << vert_indices_for_patch[3] << std::endl;
            std::vector<IndexType> reordered_vert_indices = ReorderIndicesToStartAtValue(vert_index, vert_indices_for_patch);
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: reordered_vert_indices") << reordered_vert_indices << std::endl;

            for (int i = 0; i < 4; ++i) {
                SecondRingVertices[patch_indices[i]] = reordered_vert_indices[i];
                // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: added reordered index ") << reordered_vert_indices[i] << std::endl;
            }
            // KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: SecondRingVertices ") << SecondRingVertices << std::endl;
        }
    }
    KRATOS_INFO("SecondRingVertices") << SecondRingVertices << std::endl;
    KRATOS_INFO("CATMULL_CLARK :: GetSecondRingVerticesAndFaces :: Reached end of function") << std::endl;

    return is_irregular_patch;
}

bool IsProjectionInsideCondition(NodeType& rPoint, const array_1d<double, 3>& rUnitNormal, ConditionType& rCondition, std::map<IndexType, double>& ErrorMap, Point& rLocalCoords, Point& rIntersectionPoint)
{
    // KRATOS_INFO("CATMULL_CLARK :: IsProjectionInsideCondition") << std::endl;

    // get plane description from rCondition
    // we need normal vector of the face and a point on the plane
    GeometryType& r_geom = rCondition.GetGeometry();
    
    array_1d<double, 3>& p1 = rPoint.Coordinates();
    array_1d<double, 3>& p2 = r_geom[0].Coordinates();    // coordinates of any of the nodes in condition
    array_1d<double, 3> vector_to_external_point = p1 - p2;

    // get large enough length to go from FE node to other side of plane spanned by condition geometry
    double length = 10 * norm_2(vector_to_external_point);

    array_1d<double, 3> line_start_point = rPoint.Coordinates();  // start point
    array_1d<double, 3> line_end_point;
    // double scalar_prod = inner_prod(vector_to_external_point, rUnitNormal);
    double scalar_prod = inner_prod(vector_to_external_point, r_geom.UnitNormal(0));
    if (scalar_prod > 0) {
        line_end_point = rPoint.Coordinates() - rUnitNormal * length;
    }
    else {
        line_end_point = rPoint.Coordinates() + rUnitNormal * length;
    }
    // KRATOS_INFO("IsProjectionInsideCondition :: rUnitNormal") << rUnitNormal << std::endl;
    // KRATOS_INFO("IsProjectionInsideCondition :: r_geom.UnitNormal(0)") << r_geom.UnitNormal(0) << std::endl;

    
    array_1d<double, 3> intersection_point;
    double tolerance = 0.05;

    int has_intersection = IntersectionUtilities::ComputePlaneLineIntersection(
        r_geom.Center().Coordinates(),
        r_geom.UnitNormal(0),
        line_start_point,
        line_end_point,
        intersection_point,
        tolerance
    );
    if (!has_intersection) KRATOS_ERROR << "IsProjectionInsideCondition :: No plane-line intersection found. Could not find condition in which the fe-node lies in." << std::endl;
    
    // Point intersection_point_coords;
    rIntersectionPoint.Coordinates() = intersection_point;

    // tolerance here is crucial to determine how exact the is_inside condition should be
    // the tolerance depends on how exact the fe-mesh lies on the limit
    // double tolerance = 0.05;    // tolerance for the local coordinates
    bool is_inside = r_geom.IsInside(rIntersectionPoint, rLocalCoords, tolerance);
    // KRATOS_INFO("CATMULL_CLARK :: IsProjectionInsideCondition :: tolerance ") << tolerance << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: IsProjectionInsideCondition :: rLocalCoords ") << rLocalCoords << std::endl;
    if (!is_inside) {
        IndexType current_id = rCondition.GetId();
        ErrorMap[current_id] = 0.0;

        for (auto coord_it = rLocalCoords.begin(); coord_it != rLocalCoords.end(); ++coord_it) {
            double coord = *coord_it;
            double error = abs(coord) - (1+tolerance);
            if (error > 0.0) ErrorMap[current_id] += error;
        }
    }

    return is_inside;
}

std::vector<fcpw::Vector3i> TriangulateMesh(ModelPart::MeshType& rInputMesh, std::map<IndexType, IndexType>& MapTriaIndicesToConditions, std::map<IndexType, IndexType>& MapNodeIdsToIndices)
{
    // KRATOS_INFO("CATMULL_CLARK :: FindNearestConditionOnMesh :: TriangulateMesh") << std::endl;

    std::vector<fcpw::Vector3i> triangle_indices;
    SizeType tria_index = 0;
    for (auto cond_it = rInputMesh.ConditionsBegin(); cond_it != rInputMesh.ConditionsEnd(); ++cond_it) {
        ConditionType& cond = *cond_it;
        GeometryType& geom = cond.GetGeometry();
        if (geom.size() == 3) {
            fcpw::Vector3i indices;
            indices << geom[0].GetId() - 1, geom[1].GetId() - 1, geom[2].GetId() - 1;
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

IndexType FindNearestConditionOnMesh(NodeType& rFeNode, ModelPart::MeshType& rSearchMesh)
{
    KRATOS_INFO("CATMULL_CLARK :: FindNearestConditionOnMesh") << std::endl;

    // initialize a 3d scene
    fcpw::Scene<3> scene;

    std::vector<fcpw::Vector3> mesh_nodes;
    IndexType node_index = 0;
    std::map<IndexType, IndexType> map_node_id_to_index;
    for (auto node_it = rSearchMesh.NodesBegin(); node_it != rSearchMesh.NodesEnd(); ++node_it) {
        NodeType& node = *node_it;
        fcpw::Vector3 new_node;
        new_node << node.Coordinates()[0], node.Coordinates()[1], node.Coordinates()[2];
        mesh_nodes.push_back(new_node);
        map_node_id_to_index[node.GetId()] = node_index;
        node_index += 1;
    }

    std::map<IndexType, IndexType> map_tria_index_to_condition_id;
    const std::vector<fcpw::Vector3i>& triangle_indices = TriangulateMesh(rSearchMesh, map_tria_index_to_condition_id, map_node_id_to_index);

    // load positions and indices of a single triangle mesh
    scene.setObjectCount(1);
    scene.setObjectVertices(mesh_nodes, 0);
    scene.setObjectTriangles(triangle_indices, 0);

    // build acceleration structure
    fcpw::AggregateType aggregateType = fcpw::AggregateType::Bvh_SurfaceArea;
    bool buildVectorizedBvh = true;
    scene.build(aggregateType, buildVectorizedBvh);

    // perform a closest point query
    fcpw::Interaction<3> interaction;
    fcpw::Vector3 queryPoint;
    queryPoint << rFeNode.Coordinates()[0], rFeNode.Coordinates()[1], rFeNode.Coordinates()[2];
    bool found = scene.findClosestPoint(queryPoint, interaction);

    // access distance and closest point via interaction.d and interaction.p (resp.)
    if (found) {
        int fcpw_tria_index = interaction.primitiveIndex;
        KRATOS_INFO("CATMULL_CLARK :: FindNearestConditionOnMesh :: Nearest Condition is ") << map_tria_index_to_condition_id[fcpw_tria_index] << std::endl;
        return map_tria_index_to_condition_id[fcpw_tria_index];
    }
    else KRATOS_ERROR << "FindNearestConditionOnMesh :: Could not determine nearest point on tria surface!" << std::endl;
}

IndexType NodeLiesInFace(NodeType& rFeNode, ModelPart::MeshType& rSearchMesh, Point& rIntersectionPoint)
{
    KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace") << std::endl;
    Point node_coords;
    node_coords.Coordinates() = rFeNode.Coordinates();
    // KRATOS_INFO("NodeLiesInFace :: rFeNode.Coordinates()") << rFeNode.Coordinates() << std::endl;
    // KRATOS_INFO("NodeLiesInFace :: rSearchMesh has number of conds ") << rSearchMesh.NumberOfConditions() << std::endl;

    for (auto cond_it = rSearchMesh.ConditionsBegin(); cond_it != rSearchMesh.ConditionsEnd(); ++cond_it) {
        ConditionType& r_condition = *cond_it;
        GeometryType& r_geometry = r_condition.GetGeometry();
        Point local_coords;
        double tolerance = 1e-3;
        bool is_inside = r_geometry.IsInside(node_coords, local_coords, tolerance);
        if (is_inside) {
            CoordinatesArrayType difference = node_coords.Coordinates() - r_geometry.Center();
            double diff_norm = difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2];
            if (diff_norm > 0.5) {
                // KRATOS_INFO("NodeLiesInFace :: rSearchMesh Geometry i Center ") << r_geometry.Center() << std::endl;
                // KRATOS_INFO("local_coords ") << local_coords << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Kratos condition ID ") << r_condition.GetId() << std::endl;
                continue;
            }
            return r_condition.GetId();
        }
        else {
            // KRATOS_INFO("NodeLiesInFace :: rFeNode.Coordinates()") << rFeNode.Coordinates() << std::endl;
            // KRATOS_INFO("NodeLiesInFace :: rSearchMesh Geometry i ") << r_geometry << std::endl;
            // KRATOS_INFO("local_coords ") << local_coords << std::endl;
        }
    }
    // KRATOS_ERROR;

    KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: extending search bcs FENode is not exactly on limit...") << std::endl;

    // Before failing we try to find out in which face the point should lie in.
    // Herefore we do the following:
    // 1. Search nearest neighbour (nn) in control grid with kd tree search
    // 2. Project FE node with normal of nn onto the planes of the faces adjacent to the nn
    // 3. Check all the faces if intersecting point is inside the face.
    // This should only hold true for at most one of the faces.
    // 
    // initialize nearest node search
    NodeVector control_nodes_on_limit;
    control_nodes_on_limit.reserve(rSearchMesh.Nodes().size());
    for (auto node_it = rSearchMesh.NodesBegin(); node_it != rSearchMesh.NodesEnd(); ++node_it)
    {
        control_nodes_on_limit.push_back(*(node_it.base()));
    }
    const SizeType bucket_size = 100;
    KDTree search_tree(control_nodes_on_limit.begin(), control_nodes_on_limit.end(), bucket_size);
    double distance;
    NodeTypePointer p_nearest_control_node = search_tree.SearchNearestPoint(rFeNode, distance);
    // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace : nearest control node is ") << p_nearest_control_node->GetId() << std::endl;
    // KRATOS_INFO("                                  with distance ") << distance << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace : nearest control node coords ") << p_nearest_control_node->Coordinates() << std::endl;

    std::map<IndexType, std::map<IndexType, IndexType>> map_vert_to_node, map_face_to_cond;
    OpenSubdiv::Far::TopologyRefiner * refiner = GenerateOpenSubdivRefiner(rSearchMesh, map_vert_to_node, map_face_to_cond, false);
    OpenSubdiv::Far::TopologyLevel const & firstLevel = refiner->GetLevel(0);
    IndexType nearest_control_node_osd_index = p_nearest_control_node->GetId() - 1;
    ConstIndexArray vertex_faces = firstLevel.GetVertexFaces(nearest_control_node_osd_index);
    // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace : nearest control vertex ") << nearest_control_node_osd_index << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace : GetVertexFaces(nearest_control_vertex) ") << vertex_faces[0] << " , " << vertex_faces[1] << " , " << vertex_faces[2] << " , " << vertex_faces[3] << " , " << std::endl;

    array_1d<double, 3> r_node_normal;
    for (auto face_it = vertex_faces.begin(); face_it != vertex_faces.end(); ++face_it) {
        IndexType osd_face_index = *face_it;
        IndexType kratos_cond_index = osd_face_index + 1;
        ConditionType& r_cond = rSearchMesh.GetCondition(kratos_cond_index);
        const array_1d<double, 3>& r_cond_normal_i = r_cond.GetGeometry().Normal(0);
        r_node_normal += r_cond_normal_i * r_cond.GetGeometry().Area();
    }
    r_node_normal /= norm_2(r_node_normal);
    // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace : nearest_control_vertex normal ") << r_node_normal << std::endl;
    
    std::map<IndexType, double> error_map;
    for (auto face_it = vertex_faces.begin(); face_it != vertex_faces.end(); ++face_it) {
        IndexType osd_face_index = *face_it;
        IndexType kratos_cond_index = osd_face_index + 1;
        ConditionType& r_cond = rSearchMesh.GetCondition(kratos_cond_index);
        Point local_coord;
        
        bool is_inside = IsProjectionInsideCondition(rFeNode, r_node_normal, r_cond, error_map, local_coord, rIntersectionPoint);
    
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Is projection inside condition ID ") << r_cond.GetId() << " ==> " << is_inside << std::endl;
        // KRATOS_INFO("(1 = true), (0 = false)") << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: with local coords ") << local_coord.Coordinates() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: with intersection point ") << rIntersectionPoint.Coordinates() << std::endl;

        ConstIndexArray face_verts = firstLevel.GetFaceVertices(osd_face_index);
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: GetFaceVertices(") << osd_face_index << ") : " << face_verts[0] << " , " << face_verts[1] << " , " << face_verts[2] << " , " << face_verts[3] << " , " << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Face coords ") << face_verts[0] << rSearchMesh.GetNode(face_verts[0]+1).Coordinates() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Face coords ") << face_verts[1] << rSearchMesh.GetNode(face_verts[1]+1).Coordinates() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Face coords ") << face_verts[2] << rSearchMesh.GetNode(face_verts[2]+1).Coordinates() << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: Face coords ") << face_verts[3] << rSearchMesh.GetNode(face_verts[3]+1).Coordinates() << std::endl;


        if (is_inside) {
            // KRATOS_INFO("CATMULL_CLARK :: NodeLiesInFace :: found face with nearest neighbour search for fe node id ") << rFeNode.GetId() << std::endl;
            return r_cond.GetId();
        }
    }
    // could not project the point inside a condition (face)
    // instead take the face where the projection was nearest
    IndexType index_of_smallest_error;
    double smallest_error = 100.0;
    for (auto error_it : error_map) {
        double current_error = error_it.second;
        if (current_error < smallest_error) index_of_smallest_error = error_it.first;
    }
    KRATOS_WARNING("CATMULL_CLARK :: NodeLiesInFace :: Returning the condition that is most likely to be the face the fe node lies in") << smallest_error << std::endl;
    return index_of_smallest_error;

    // KRATOS_ERROR << "NodeLiesInFace :: Could not find condition in which the fe-node lies in. \nHint: Check if the FE mesh lies on the limit surface." << std::endl;
}

std::vector<double> GetLimitWeights(SizeType NeighbourhoodSize)
{
    // The weights for computing the limit stem from an analytical approach 
    // where the spectral analysis of the local subdivision matrix is used to 
    // determine the weights in the limit of every point in the neighbourhood.
    // This enables to compute every point on the limit, even irregular points.
    if (NeighbourhoodSize == 7)
        {
            return {0.375, 1.0/6.0, 1.0/24.0, 1.0/6.0, 1.0/24.0, 1.0/6.0, 1.0/24.0};
        }
        else if (NeighbourhoodSize == 9)    // regular point for Catmull-Clark
        {
            return {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};
        }
        else if (NeighbourhoodSize == 11)
        {
            return {0.5, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02, 0.08, 0.02};
        }
        else if (NeighbourhoodSize == 12)
        {
            return {6.0/11.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 
            5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0, 2.0/33.0, 5.0/333.0};
        }
        else {
            std::vector<double> vec;
            vec.clear();
            return vec;
        }
}

void PushNodesToLimit(ModelPart& rInputModelPart, ModelPart& rOutputModelPart, std::map<IndexType, std::vector<IndexType>> rFirstRingNodeNeighbourhoods)
{   
    KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit") << std::endl;
    // #pragma omp parallel
    // for (SizeType i = 0; i < rInputModelPart.NumberOfNodes(); ++i)
    for (auto node_it = rInputModelPart.NodesBegin(); node_it != rInputModelPart.NodesEnd(); ++node_it)
    {
        const NodeType& r_node_i = *node_it;
        IndexType node_id = r_node_i.GetId();
        // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit :: Current Kratos node ID") << node_id << std::endl;
        // const CoordinatesArrayType& r_coords = r_node_i.Coordinates();

        std::vector<IndexType> r_node_i_neighbour_indices = rFirstRingNodeNeighbourhoods[node_id];
        // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit :: with first ring neighbourhood") << r_node_i_neighbour_indices << std::endl;
        
        SizeType neighb_size = r_node_i_neighbour_indices.size();
        std::vector<CoordinatesArrayType> neighbourhood_coords;
        
        for (auto neighb_idx_it = r_node_i_neighbour_indices.begin(); neighb_idx_it != r_node_i_neighbour_indices.end(); ++neighb_idx_it)
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
        if (limit_weights.size() == 0) {
            rOutputModelPart.RemoveNode(node_id);
        }
        else {
            rOutputModelPart.GetNode(node_id).Coordinates()[0] = coords[0];
            rOutputModelPart.GetNode(node_id).Coordinates()[1] = coords[1];
            rOutputModelPart.GetNode(node_id).Coordinates()[2] = coords[2];
        }
    }
    // KRATOS_INFO("CATMULL_CLARK :: PushNodesToLimit successful") << std::endl;
}

std::vector<double> GetUnitWeightDistributionForNode(IndexType osd_v_idx, OpenSubdiv::Far::TopologyRefiner * refiner)
{
    // KRATOS_INFO("CATMULL_CLARK :: GetUnitWeightDistributionForNode") << std::endl;

    // Create a buffer to hold the position of the refined verts and
    // local points, then copy the coarse positions at the beginning.
    
    SizeType n_uniformRef = refiner->GetNumLevels() - 1;
    SizeType nverts = refiner->GetLevel(0).GetNumVertices();

    std::vector<double> vertices(3*nverts, 0.0);
    std::vector<double> cp_weights(nverts, 0.0);
    cp_weights[osd_v_idx] = 1.0;

    Vertex *verts = (Vertex*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(Vertex) );
    memcpy(&verts[0], &vertices[0], nverts*sizeof(double));

    VertexValue *vverts = (VertexValue*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(VertexValue) );
    memcpy(&vverts[0], &cp_weights[0], nverts*sizeof(double));

    // Interpolate vertex primvar data at each level
    OpenSubdiv::Far::PrimvarRefiner primvarRefiner(*refiner);

    Vertex * src = verts;
    VertexValue * vsrc = vverts;
    for (SizeType level = 1; level <= n_uniformRef; ++level) {
        Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;

        VertexValue * vdst = vsrc + refiner->GetLevel(level-1).GetNumVertices();
        primvarRefiner.Interpolate(level, vsrc, vdst);
        vsrc = vdst;
    }

    // return the current vertices and faces
    OpenSubdiv::Far::TopologyLevel const & refLastLevel = refiner->GetLevel(n_uniformRef);
    
    nverts = refLastLevel.GetNumVertices();
    
    // copy vertex information
    int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
    vertices.resize(nverts*3);
    memcpy(&vertices[0], &verts[firstOfLastVerts], nverts*3*sizeof(double));

    cp_weights.resize(nverts);
    memcpy(&cp_weights[0], &vverts[firstOfLastVerts], nverts*sizeof(double));

    // KRATOS_INFO("CATMULL_CLARK :: GetUnitWeightDistributionForNode :: cp_weights") << cp_weights << std::endl;

    return cp_weights;
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
        "entity_type"                                 : "condition",
        "output_path"                                 : "VTK_Output"
    })" );
    VtkOutput(rOuptutMP, parameters).PrintOutput(filename);
}

std::vector<NodeTypePointer> CreateRegularPatchForFeNode(
    NodeType& rFeNode, 
    IndexType KratosFaceId, 
    ModelPart::MeshType& rControlGrid,
    ModelPart::MeshType& rInputMesh, 
    ModelPart& rInputModelPart, 
    std::map<IndexType, Vector>& RefiningWeights,
    std::map<IndexType, Point>& rLimitIntersectionPointsMap,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_vert_to_node,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_face_to_cond,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode in Kratos Face ") << KratosFaceId << std::endl;
    auto geom = rInputMesh.GetCondition(KratosFaceId).GetGeometry();
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom IDs ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;

    std::vector<IndexType> second_ring_vertices;    // contains osd ids
    std::vector<IndexType> second_ring_faces;       // contains osd ids
    bool patch_is_irregular = GetSecondRingVerticesAndFaces(KratosFaceId, rInputMesh, second_ring_vertices, second_ring_faces);
    KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: second_ring_vertices ") << second_ring_vertices << std::endl;
    KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: second_ring_faces ") << second_ring_faces << std::endl;

    SizeType max_refining_depth = 4;
    SizeType current_refining_depth = 1;    // we start at first refinement of control polygon
    Model patch_model;
    Model refined_patch_model;
    Model refined_limit_model;
    // ModelPart& r_patch_mp = refined_patch_model.CreateModelPart("patch", 3);
    ModelPart& r_refined_patch_mp = refined_patch_model.CreateModelPart("refined_patch", 3);
    ModelPart& r_refined_limit_mp = refined_limit_model.CreateModelPart("refined_limit", 3);
    // r_patch_mp.CreateNewProperties(0);
    r_refined_patch_mp.CreateNewProperties(0);
    r_refined_limit_mp.CreateNewProperties(0);
    std::map<IndexType, std::vector<IndexType>> first_ring_vert_map;

    // r_patch_mp.SetNodes(rInputMesh.pNodes());
    // r_patch_mp.SetConditions(rInputMesh.pConditions());
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;
    // ModelPart& r_patch_mp = AuxiliarModelPartUtilities(rInputModelPart).DeepCopyModelPart("patch", &patch_model);
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: after DeepCopyModelPart of modelpart with name ") << rInputModelPart.Name() << std::endl;

    // NodeTypePointer as return type
    std::vector<NodeTypePointer> rPatchNodesVector;
    // refine the patch until fe-node is inside regular patch
    while (patch_is_irregular && max_refining_depth > current_refining_depth) 
    {
        current_refining_depth += 1;
        KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: current_refining_depth") << current_refining_depth << std::endl;

        // auto geom = rInputMesh.GetCondition(KratosFaceId).GetGeometry();
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom ") << geom << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: geom IDs ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;

        // CreateModelPartFor2RingPatch(rInputMesh, r_patch_mp, second_ring_vertices, second_ring_faces);   // for now we refine full geometry to get local refinement, can be improved by working on 2-ring patch only
        
        RefineGeometry(rInputModelPart.GetMesh(), r_refined_patch_mp, r_global_map_vert_to_node, r_global_map_face_to_cond, false, current_refining_depth, refiner);
        // RefineGeometry(r_patch_mp.GetMesh(), r_refined_patch_mp, false, current_refining_depth, refiner);
        CreateFirstRingNeighbourhoodsOfNodes(r_refined_patch_mp.GetMesh(), first_ring_vert_map, r_global_map_vert_to_node, r_global_map_face_to_cond, refiner);
        // create modelpart whose nodes lie on the limit surface
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: DeepCopyModelPart of modelpart with name ") << r_refined_patch_mp.Name() << std::endl;
        ModelPart& r_refined_limit_mp = AuxiliarModelPartUtilities(r_refined_patch_mp).DeepCopyModelPart(r_refined_patch_mp.Name()+"_limit_"+std::to_string(current_refining_depth), &refined_limit_model);
        PushNodesToLimit(r_refined_patch_mp, r_refined_limit_mp, first_ring_vert_map);
        
        // after refining: check which face fe_node lies in
        // initialize nearest node search
        // NodeVector control_nodes_on_limit;
        // control_nodes_on_limit.reserve(r_refined_limit_mp.Nodes().size());
        // for (auto node_it = r_refined_limit_mp.NodesBegin(); node_it != r_refined_limit_mp.NodesEnd(); ++node_it)
        // {
        //     control_nodes_on_limit.push_back(*(node_it.base()));
        // }
        // const SizeType bucket_size = 100;
        // KDTree search_tree(control_nodes_on_limit.begin(), control_nodes_on_limit.end(), bucket_size);

        // get nearest node
        // double distance;
        // currently not necessary
        // NodeTypePointer p_nearest_control_node = search_tree.SearchNearestPoint(rFeNode, distance);
        // identify the face that the node lies in
        // Point intersection_point;
        // KratosFaceId = NodeLiesInFace(rFeNode, r_refined_limit_mp.GetMesh(), intersection_point);
        KratosFaceId = FindNearestConditionOnMesh(rFeNode, r_refined_limit_mp.GetMesh());
        // rLimitIntersectionPointsMap[rFeNode.GetId()] = intersection_point;
        // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New KratosFaceId in refined mesh ") << KratosFaceId << std::endl;
        std::map<IndexType, IndexType> r_map_face_to_cond = r_global_map_face_to_cond[current_refining_depth];

        // create 2-ring vertex neighbourhoods around the face the fe-node lies in
        patch_is_irregular = GetSecondRingVerticesAndFaces(KratosFaceId, r_refined_patch_mp.GetMesh(), second_ring_vertices, second_ring_faces);
        KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: New patch for KratosFaceId: patch_is_irregular ") << patch_is_irregular << std::endl;
        if (!patch_is_irregular) {
            std::string filename = r_refined_patch_mp.Name() + std::to_string(rFeNode.GetId());
            OutputModelPartToVtk(r_refined_patch_mp, filename);
            // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: refined regular patch (osd ids) is ") << second_ring_vertices << std::endl;
            for (IndexType i = 0; i < second_ring_vertices.size(); ++i) {
                IndexType kratos_id = second_ring_vertices[i] + 1;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode osd_id ") << second_ring_vertices[i] << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode kratos_id ") << kratos_id << std::endl;
                // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode r_refined_patch_mp.GetNode(kratos_id) ") << r_refined_patch_mp.GetNode(kratos_id) << std::endl;
                rPatchNodesVector.push_back(r_refined_patch_mp.pGetNode(kratos_id));
            }

            // get weights from refining: meaning the influence of the base control points on the refined vertices of the regular patch
            IndexType osd_face_index = GetOSDIndexFromKratosID(KratosFaceId, r_map_face_to_cond);

            // if (!refiner) {
            //     refiner = GenerateOpenSubdivRefiner(r_patch_mp.GetMesh(), false);
            // }
            // refiner->RefineUniform(OpenSubdiv::Far::TopologyRefiner::UniformOptions(current_refining_depth));
            // TODO: only refiner if level does not exost yet, avoids to refine every time anew
            OpenSubdiv::Far::TopologyLevel const & lastRefLevel = refiner->GetLevel(current_refining_depth);
            // get face index of control grid where FE node lies in
            for (SizeType ref = 0; ref < current_refining_depth; ++ref) {
                osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
            }
            // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: Parent osd_face_index on control polygon ") << osd_face_index << std::endl;
            IndexType kratos_base_face_index = r_map_face_to_cond[osd_face_index];

            std::vector<IndexType> base_second_ring_vertices;    // contains osd ids
            std::vector<IndexType> base_second_ring_faces;       // contains osd ids

            // get second ring verts and faces for the base control grid P0
            GetSecondRingVerticesAndFaces(kratos_base_face_index, rControlGrid, base_second_ring_vertices, base_second_ring_faces);
            // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: base_second_ring_vertices ") << base_second_ring_vertices << std::endl;

            // now base_second_ring_vertices contains all osd ids of the ctrl pts that influence the vertices of final regular patch
            // evaluate unit weight distribution for every node on the base level (refinement 0)

            for(auto v_it = base_second_ring_vertices.begin(); v_it != base_second_ring_vertices.end(); ++v_it) {
                Vector influence_weights_of_base_on_refined(16);
                IndexType osd_v_id = *v_it;
                IndexType kratos_v_id = osd_v_id + 1;
                
                std::vector<double> weight_distr_of_cp_i = GetUnitWeightDistributionForNode(osd_v_id, refiner);

                IndexType idx = 0;
                for(auto regular_patch_it = second_ring_vertices.begin(); regular_patch_it != second_ring_vertices.end(); ++ regular_patch_it) {
                    IndexType osd_patch_vertex_id = *regular_patch_it;
                    influence_weights_of_base_on_refined[idx] = weight_distr_of_cp_i[osd_patch_vertex_id];
                    idx += 1;
                }

                RefiningWeights[kratos_v_id] = influence_weights_of_base_on_refined;
            }
        }
    }
    if (current_refining_depth >= max_refining_depth && patch_is_irregular) {
        // throw error
        KRATOS_ERROR << "Reached maximum refining depth, could not determine regular patch for FE-node " << rFeNode.GetId() << std::endl;
    }

    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: rPatchNodesVector ") << rPatchNodesVector << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateRegularPatchForFeNode :: end of function ") << std::endl;
    return rPatchNodesVector;
}

void CreateFirstRingNeighbourhoodsOfConditions(
    ModelPart::MeshType& rInputMesh,
    std::map<IndexType, std::vector<IndexType>>& rOutputMap,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_verts,
    std::map<IndexType, std::map<IndexType, IndexType>>& r_global_map_faces,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    /*
    Structuring/numbering of the neighbourhood's nodes should be according to NURBS surface evaluation

        1----->2----->3----->4
        |      |      |      |
        |      |      |      |
        5----->6----->7----->8
        |      |      |      |
        |      |      |      |
        9----->10---->11---->12
        |      |      |      |
        |      |      |      |
        13---->14---->15---->16

    The regular patches will only be created for quad elements with regular vertices.
    */
    KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfConditions ") << std::endl;

    SizeType ref_level = 0;

    // if no refiner passed, generate new one based on rInputMesh
    if (!refiner) {
        /// generate opensubdiv refiner
        refiner = GenerateOpenSubdivRefiner(rInputMesh, r_global_map_verts, r_global_map_faces, false);
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfConditions :: creating new refiner ") << std::endl;
    }
    else {
        ref_level = refiner->GetNumLevels() - 1;
        // KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfConditions :: using specified refiner ") << std::endl;
    }
    // get topology from level ref_level
    OpenSubdiv::Far::TopologyLevel const & refLevel = refiner->GetLevel(ref_level);
    SizeType nfaces = refLevel.GetNumFaces();

    std::vector<std::vector<IndexType>> patchPointsIndicesPerFace = { {  5,  4,  0,  1 },
                                                                      {  6,  2,  3,  7 },
                                                                      { 10, 11, 15, 14 },
                                                                      {  9, 13, 12,  8 } };

    int regular_patches = 0;
    int irregular_patches = 0;

    std::map<IndexType, IndexType> map_verts = r_global_map_verts[ref_level];
    std::map<IndexType, IndexType> map_faces = r_global_map_faces[ref_level];

    for (IndexType face_i = 0; face_i < nfaces; ++face_i) {
        std::vector<IndexType> regular_patch_indices(16);

        IndexType condition_id = map_faces[face_i];
        bool irregular_patch = false;
        
        SizeType num_verts = refLevel.GetFaceVertices(face_i).size();
        if (num_verts != 4) irregular_patch = true;

        ConstIndexArray vert_indices = refLevel.GetFaceVertices(face_i);
        for (auto vert_index_it = vert_indices.begin(); vert_index_it != vert_indices.end(); ++vert_index_it) {
            IndexType vert_index = *vert_index_it;
            if(!refLevel.IsVertexValenceRegular(vert_index)) irregular_patch = true;

            ConstIndexArray face_indices = refLevel.GetVertexFaces(vert_index);
            for (auto face_index_it = face_indices.begin(); face_index_it != face_indices.end(); ++face_index_it) {
                IndexType face_index = *face_index_it;
                num_verts = refLevel.GetFaceVertices(face_index).size();
                if (num_verts != 4) irregular_patch = true;
            }
        }

        if(irregular_patch) {
            rOutputMap[condition_id] = {};
            irregular_patches += 1;
            continue;
        }

        IndexType patch_position_index = 0;
        for (auto vert_index_it = vert_indices.begin(); vert_index_it != vert_indices.end(); ++vert_index_it) {
            IndexType vert_index = *vert_index_it;
            std::vector<IndexType> patch_indices = patchPointsIndicesPerFace[patch_position_index];
            patch_position_index += 1;

            ConstIndexArray face_indices = refLevel.GetVertexFaces(vert_index);
            IndexType opposite_face_index = GetOppositeIndex(face_indices, face_i);
            ConstIndexArray vert_indices_for_patch = refLevel.GetFaceVertices(opposite_face_index);
            std::vector<IndexType> reordered_vert_indices = ReorderIndicesToStartAtValue(vert_index, vert_indices_for_patch);

            for (int i = 0; i < 4; ++i) {
                regular_patch_indices[patch_indices[i]] = reordered_vert_indices[i];
            }
        }
        rOutputMap[condition_id] = regular_patch_indices;
        regular_patches += 1;
    }
    KRATOS_INFO("CATMULL_CLARK :: CreateFirstRingNeighbourhoodsOfConditions reached end of function") << std::endl;
}

// not needed?
// double ComputeAreaOfTriangle(NodeTypePointer pNode0, NodeType& rNode1, NodeType& rNode2) 
// {
//     const array_1d<double, 3>& P = pNode0->Coordinates();
//     const array_1d<double, 3>& A = rNode1.Coordinates();
//     const array_1d<double, 3>& B = rNode2.Coordinates();

//     const array_1d<double, 3> PA = A - P;
//     const array_1d<double, 3> PB = B - P;

//     array_1d<double, 3> cross_product;
//     MathUtils<double>::CrossProduct(cross_product, PA, PB);
//     double cp0 = cross_product[0];
//     double cp1 = cross_product[1];
//     double cp2 = cross_product[2];
    
//     return 0.5 * sqrt(cp0*cp0 + cp1*cp1 + cp2*cp2);
// };

// not needed?
// bool CheckIfRegularPatch(std::vector<IndexType> Patch, ModelPart::MeshType& rInputMesh) 
// {
//     bool is_regular = true;
//     std::map<IndexType, std::map<IndexType, IndexType>> map_vert, map_face;
//     OpenSubdiv::Far::TopologyRefiner* refiner = GenerateOpenSubdivRefiner(rInputMesh, map_vert, map_face, false);
//     OpenSubdiv::Far::TopologyLevel const & firstLevel = refiner->GetLevel(0);
    
//     SizeType nverts = firstLevel.GetNumVertices();
//     for (IndexType i = 0; i < nverts; ++i) {
//         is_regular = firstLevel.IsVertexValenceRegular(i);
//     }

//     return is_regular;
// };

// }   // namespace SDSUtils


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

    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: original number of Control Nodes ") << rControlPolygon.NumberOfNodes() << std::endl;
    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: original number of Control Conds ") << rControlPolygon.NumberOfConditions() << std::endl;


    Model model;
    Model model_limit;
    ModelPart& r_first_subd_for_quads = model.CreateModelPart("first_subd", 3);
    r_first_subd_for_quads.CreateNewProperties(0);

    // the maps store for every refinement level the relation between Kratos IDs and OSD indices
    // structure is the following:
    // map<level, map<KratosID, OSDIndex>>
    std::map<IndexType, std::map<IndexType, IndexType>> global_map_osd_vert_to_kratos_node_ids;
    std::map<IndexType, std::map<IndexType, IndexType>> global_map_osd_face_to_kratos_cond_ids;
    /// generate opensubdiv refiner
    OpenSubdiv::Far::TopologyRefiner * refiner_main = GenerateOpenSubdivRefiner(rControlPolygon.GetMesh(), global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, FixFreeEdges);

    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: refiner_main->GetNumVerticesTotal() ") << refiner_main->GetNumVerticesTotal() << std::endl;
    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: refiner_main->GetNumFacesTotal() ") << refiner_main->GetNumFacesTotal() << std::endl;

    RefineGeometry(rControlPolygon.GetMesh(), r_first_subd_for_quads, global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, FixFreeEdges, 1, refiner_main);

    // create first ring neighbourhood of nodes of control polygon
    CreateFirstRingNeighbourhoodsOfNodes(r_first_subd_for_quads.GetMesh(), mFirstRingNodes, global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, refiner_main);  // enables pushing the center node to the limit
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: after mFirstRingNodes.size ") << mFirstRingNodes.size() << std::endl;
    // create first ring neighbourhood of elements of control polygon
    CreateFirstRingNeighbourhoodsOfConditions(r_first_subd_for_quads.GetMesh(), mFirstRingFaces, global_map_osd_vert_to_kratos_node_ids, global_map_osd_face_to_kratos_cond_ids, refiner_main);
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: mFirstRingFaces.size ") << mFirstRingFaces.size() << std::endl;
    
    // create modelpart whose nodes lie on the limit surface
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: DeepCopyModelPart of modelpart with name ") << r_first_subd_for_quads.Name() << std::endl;
    ModelPart& r_points_on_limit = AuxiliarModelPartUtilities(r_first_subd_for_quads).DeepCopyModelPart(r_first_subd_for_quads.Name()+"_limit", &model_limit);

    PushNodesToLimit(r_first_subd_for_quads, r_points_on_limit, mFirstRingNodes);

    /// for debugging purposes start //
    std::map<IndexType, std::vector<IndexType>> base_cp_first_ring_neighb;
    std::map<IndexType, std::map<IndexType, IndexType>> vert_map, face_map;
    CreateFirstRingNeighbourhoodsOfConditions(rControlPolygon.GetMesh(), base_cp_first_ring_neighb, vert_map, face_map);

    // loop over all nodes and check their read ids
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rControlPolygon Node IDs : ") << std::endl;
    // IndexType index = 0;
    // for(auto node_it = rControlPolygon.NodesBegin(); node_it != rControlPolygon.NodesEnd(); ++node_it) {
    //     NodeType& r_node = *node_it;
    //     IndexType r_node_id = r_node.GetId();
    //     KRATOS_INFO("#") << index << " ~ " << r_node_id << " ";
    //     index += 1;
    // }
    // KRATOS_INFO(".") << std::endl;
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rControlPolygon Node IDs end") << std::endl;
    // KRATOS_ERROR << "Error out for debug" << std::endl;
    /// for debugging purposes end //


    // // initialize nearest node search
    // NodeVector control_nodes_on_limit;
    // control_nodes_on_limit.reserve(r_points_on_limit.Nodes().size());
    // for (auto node_it = r_points_on_limit.NodesBegin(); node_it != r_points_on_limit.NodesEnd(); ++node_it)
    // {
    //     control_nodes_on_limit.push_back(*(node_it.base()));
    // }
    // const SizeType bucket_size = 100;
    // KDTree search_tree(control_nodes_on_limit.begin(), control_nodes_on_limit.end(), bucket_size);
    
    // mapping matrix from control points to fe points is mMappingMatrix
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: check before loop ") << std::endl;

    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: total number of FE Nodes ") << rControlledMesh.NumberOfNodes() << std::endl;

    // save local coordinates of fe nodes that are not exactly on limit
    // maps an fe node's id to the corresponding local coordinates on the limit
    std::map<IndexType, Point> limit_intersection_points_map;

    // for (SizeType i = 0; i < rControlledMesh.NumberOfNodes(); i++) {
    IndexType fe_node_index = 0;
    for (auto fe_node_it = rControlledMesh.Nodes().begin(); fe_node_it != rControlledMesh.Nodes().end(); ++fe_node_it) {
        IndexType row_index = fe_node_index * rControlPolygon.NumberOfNodes();
        fe_node_index += 1;
        NodeType& fe_node = *fe_node_it;
        IndexType fe_node_id = fe_node.GetId();
        KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: starting with fe node id ") << fe_node_id << std::endl;

        std::vector<IndexType> patch_ids;
        std::vector<NodeTypePointer> patch;
        // 1. check if fe point lies in a regular patch for parameterizing the limit surface it lies on
        //      1.1 create neighbourhood for regular patch in the 2nd neighbouring ring 
        //          start by refining two times to have at most 1 irregular point in a regular patch
        //          1.1.1 push all control nodes to limit surface

        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: check inside loop 1 ") << std::endl;
        //          1.1.2 first find nearest neighbour P1 on fe mesh then find element/face the point lies in                    
                    // nearest neighbour search currently not necessary since all faces are checked if fe-node is inside
                    // double distance;
                    // NodeTypePointer p_nearest_control_node = search_tree.SearchNearestPoint(fe_node, distance);
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: distance ") << distance << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: check inside loop 2 ") << std::endl;
                    // check which face the fe_node lies in for regular patch creation
                    // Point intersection_point;
                    // IndexType kratos_face_id = NodeLiesInFace(fe_node, r_points_on_limit.GetMesh(), intersection_point);
                    IndexType kratos_face_id = FindNearestConditionOnMesh(fe_node, r_points_on_limit.GetMesh());
                    // limit_intersection_points_map[fe_node.GetId()] = intersection_point;

        //          1.1.3 find neighbours of P1 - P2, P3, ... and create regular patch neighbourhood
                    patch_ids = mFirstRingFaces[kratos_face_id];
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: patch_ids.size() ") << patch_ids.size() << std::endl;
                    for (auto node_it = patch_ids.begin(); node_it != patch_ids.end(); ++node_it) {
                        IndexType osd_vertex_index = *node_it;
                        IndexType kratos_node_id = global_map_osd_vert_to_kratos_node_ids[1][osd_vertex_index];
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: *node_it ") << *node_it << std::endl;
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: kratos_node_id ") << kratos_node_id << std::endl;
                        patch.push_back(r_first_subd_for_quads.pGetNode(kratos_node_id));
                    }
        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: Before calling CreateRegularPatchForFeNode ") << std::endl;
        //          1.1.4 if any Pi has valence different to 4 then subdivide the original control neighbourhood of Pi and jump to 1.1.1 
                    bool is_irregular = (patch_ids.size() == 0);
                    // create refining weight values for refined regular patches if patch is irregular
                    std::map<IndexType, Vector> refining_weights_map;
                    if(is_irregular) {
                        
                        // create regular patch by refinement
                        patch = CreateRegularPatchForFeNode(
                            fe_node, 
                            kratos_face_id, 
                            rControlPolygon.GetMesh(), 
                            r_first_subd_for_quads.GetMesh(), 
                            rControlPolygon, 
                            refining_weights_map, 
                            limit_intersection_points_map, 
                            global_map_osd_vert_to_kratos_node_ids, 
                            global_map_osd_face_to_kratos_cond_ids,
                            refiner_main
                        );
                            // r_first_subd_for_quads, 
                        // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: returned patch ") << patch << std::endl;

                        // nur vorübergehend, um r_first_subd_for_quads zu füllen
                        // RefineGeometry(rControlPolygon.GetMesh(), r_first_subd_for_quads, FixFreeEdges);
                    }
                    else {
                        // get refining weights
                        IndexType osd_face_index = global_map_osd_vert_to_kratos_node_ids[1][kratos_face_id];
                        OpenSubdiv::Far::TopologyLevel const & lastRefLevel = refiner_main->GetLevel(1);
                        osd_face_index = lastRefLevel.GetFaceParentFace(osd_face_index);
                        IndexType kratos_base_face_index = osd_face_index + 1;
                        std::vector<IndexType> base_second_ring_vertices;    // contains osd ids
                        std::vector<IndexType> base_second_ring_faces;       // contains osd ids
                        bool test_irregular = GetSecondRingVerticesAndFaces(kratos_base_face_index, rControlPolygon.GetMesh(), base_second_ring_vertices, base_second_ring_faces);
                        
                        KRATOS_INFO("test_irregular") << test_irregular << std::endl;
                        KRATOS_INFO("base_second_ring_vertices") << base_second_ring_vertices << std::endl;
                        KRATOS_INFO("base_second_ring_faces") << base_second_ring_faces << std::endl;
                        // KRATOS_INFO("global_map_osd_vert_to_kratos_node_ids[0]") << std::endl;
                        // for (auto it = global_map_osd_vert_to_kratos_node_ids[0].begin(); it != global_map_osd_vert_to_kratos_node_ids[0].end(); ++it) {
                        //     KRATOS_INFO(" / osd vert idx") << it->first << " kratos node id: " << it->second;
                        // }

                        for(auto v_it = base_second_ring_vertices.begin(); v_it != base_second_ring_vertices.end(); ++v_it) {
                            Vector influence_weights_of_base_on_refined(16);
                            IndexType osd_v_id = *v_it;
                            IndexType kratos_v_id = global_map_osd_vert_to_kratos_node_ids[0][osd_v_id];
                            
                            std::vector<double> weight_distr_of_cp_i = GetUnitWeightDistributionForNode(osd_v_id, refiner_main);

                            IndexType idx = 0;
                            for(auto regular_patch_it = patch_ids.begin(); regular_patch_it != patch_ids.end(); ++ regular_patch_it) {
                                IndexType osd_patch_vertex_id = *regular_patch_it;
                                influence_weights_of_base_on_refined[idx] = weight_distr_of_cp_i[osd_patch_vertex_id];
                                idx += 1;
                            }

                            refining_weights_map[kratos_v_id] = influence_weights_of_base_on_refined;
                        }
                    }
        

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
                p_surface_geometry->ProjectionPointGlobalToLocalSpace(fe_node.Coordinates(), local_coords);
                // if (limit_intersection_points_map.find(fe_node.GetId()) == limit_intersection_points_map.end()) {
                //     p_surface_geometry->ProjectionPointGlobalToLocalSpace(fe_node.Coordinates(), local_coords);
                // }
                // else {
                //     CoordinatesArrayType intersection_coords = limit_intersection_points_map[fe_node.GetId()].Coordinates();
                //     p_surface_geometry->ProjectionPointGlobalToLocalSpace(intersection_coords, local_coords);
                // }
                KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: FE global_coords   ") << fe_node.Coordinates() << std::endl;
                KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: FE local_coords    ") << local_coords << std::endl;

        // 3. evaluate the fe point on limit surface to establish relation to (weights of) the control points
                Vector r_weights;
                p_surface_geometry->ShapeFunctionsValues(r_weights, local_coords);
                CoordinatesArrayType global_coords;
                global_coords[0] = 0.0;
                global_coords[1] = 0.0;
                global_coords[2] = 0.0;
                for (SizeType i = 0; i < r_weights.size(); ++i) {
                    global_coords[0] += r_weights[i] * points[i].Coordinates()[0];
                    global_coords[1] += r_weights[i] * points[i].Coordinates()[1];
                    global_coords[2] += r_weights[i] * points[i].Coordinates()[2];
                }
                KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: recomputed global_coords") << global_coords << std::endl;
                KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: control point weights ") << r_weights << std::endl;
                CoordinatesArrayType difference;
                difference = fe_node.Coordinates() - global_coords;
                double diff_norm = difference[0]*difference[0] + difference[1]*difference[1] + difference[2]*difference[2];
                diff_norm = sqrt(diff_norm);
                KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: diff_norm ") << diff_norm << std::endl;
                if (diff_norm > 0.5) {
                    KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: control points ") << points << std::endl;
                }

        //      3.1 save the relation in mapping matrix
                // todo : treat refined regular patches
                if (is_irregular) {
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: refining_weights_map.size() ") << refining_weights_map.size() << std::endl;

                    for (auto ctrl_pt_it : refining_weights_map) {
                        // KRATOS_INFO(" ") << ctrl_pt_it.first << " : " << ctrl_pt_it.second << std::endl;

                        IndexType ctrl_node_id = ctrl_pt_it.first;
                        // KRATOS_INFO(" row_index ") << row_index << std::endl;
                        // KRATOS_INFO(" control_node_id ") << ctrl_node_id << std::endl;
                        // KRATOS_INFO(" entry index = row_index + control_node_id-1 ") << row_index + ctrl_node_id-1 << std::endl;
                        // KRATOS_INFO(" inner_prod(r_weights, ctrl_pt_it.second) ") << inner_prod(r_weights, ctrl_pt_it.second) << std::endl;
                        // output should have data structure: matrix dim = [num_fe_nodes x num_ctrl_nodes]
                        // thus every row takes num_ctrl_pts entries
                        // since we do outer loop over fe nodes, we need to set entries for the column here
                        rOutputData[row_index + ctrl_node_id-1] = inner_prod(r_weights, ctrl_pt_it.second);
                    }
                }
                else {
                    for (auto ctrl_pt_it : refining_weights_map) {
                        KRATOS_INFO(" ") << ctrl_pt_it.first << " : " << ctrl_pt_it.second << std::endl;

                        IndexType ctrl_node_id = ctrl_pt_it.first;
                        KRATOS_INFO(" row_index ") << row_index << std::endl;
                        KRATOS_INFO(" control_node_id ") << ctrl_node_id << std::endl;
                        KRATOS_INFO(" entry index = row_index + control_node_id-1 ") << row_index + ctrl_node_id-1 << std::endl;
                        KRATOS_INFO(" inner_prod(r_weights, ctrl_pt_it.second) ") << inner_prod(r_weights, ctrl_pt_it.second) << std::endl;
                        // output should have data structure: matrix dim = [num_fe_nodes x num_ctrl_nodes]
                        // thus every row takes num_ctrl_pts entries
                        // since we do outer loop over fe nodes, we need to set entries for the column here
                        rOutputData[row_index + ctrl_node_id-1] = inner_prod(r_weights, ctrl_pt_it.second);
                    }
                    // for (auto node_it = patch.begin(); node_it != patch.end(); ++node_it) {
                    //     const NodeTypePointer p_node = *node_it;
                    //     SizeType number_of_control_points = rControlPolygon.NumberOfNodes();
                    //     for(IndexType w_idx = 0; w_idx < r_weights.size(); ++w_idx) {
                    //         // KRATOS_INFO(" writing into ") << r_weights << std::endl;
                    //         rOutputData[row_index + p_node->GetId()-1] = r_weights[w_idx];   // GetId might not always be the best way, since some solvers use very high ID values
                    //     }
                    // }
                }
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rOutputData ") << rOutputData << std::endl;
                
                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: fe_node_index ") << fe_node_index << std::endl;

                // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: row_index ") << row_index << std::endl;


                bool error_out_in_first_loop = false;
                KRATOS_ERROR_IF(error_out_in_first_loop);
                // KRATOS_ERROR_IF(fe_node_index == 3);

        KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: successfully finished loop number ") << fe_node.GetId() << std::endl;

    }
    // KRATOS_INFO("CATMULL_CLARK :: CreateMappingMatrix :: rOutputData ") << rOutputData << std::endl;
    KRATOS_INFO("Finished creating mapping matrix...") << std::endl;
};

} // namespace kratos