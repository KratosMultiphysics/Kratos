//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Bastian Devresse,
//

#ifndef OPENSUBDIV_UTILITIES_H
#define OPENSUBDIV_UTILITIES_H

// Project includes
#include "includes/model_part.h"
#include "includes/mesh.h"
#include "includes/condition.h"

// External includes
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::seconds;
/// OpenSubdiv standard paths
// #include <opensubdiv/far/topologyDescriptor.h>
// #include <opensubdiv/far/primvarRefiner.h>
// #include <opensubdiv/far/patchTableFactory.h>
// #include <opensubdiv/far/patchMap.h>
/// OpenSubdiv specific paths
// #include "../../../custom_external_libraries/OpenSubdiv/opensubdiv/far/topologyDescriptor.h"
// #include "../../../custom_external_libraries/OpenSubdiv/opensubdiv/far/primvarRefiner.h"
// #include "../../../custom_external_libraries/OpenSubdiv/opensubdiv/far/patchTableFactory.h"
// #include "../../../custom_external_libraries/OpenSubdiv/opensubdiv/far/patchMap.h"
/// OpenSubdiv specific paths
#include "opensubdiv/far/topologyDescriptor.h"
#include "opensubdiv/far/primvarRefiner.h"
#include "opensubdiv/far/patchTableFactory.h"
#include "opensubdiv/far/patchMap.h"

/// This file contains classes for handling OpenSubdiv's types
/// and functions for applying OpenSubdivs utilities

using namespace Kratos;

using IndexType = std::size_t;
using SizeType = std::size_t;
using ConditionType = ModelPart::ConditionType;
using GeometryType = Geometry<Node>;
using MeshType = ModelPart::MeshType;
using NodesContainerType = ModelPart::NodesContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;
using NodesArrayType = Condition::NodesArrayType;
typedef NodeType Node;
typedef NodeType::Pointer NodeTypePointer;
typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;
typedef OpenSubdiv::Far::TopologyDescriptor Descriptor;

// Vertex container implementation.
struct Vertex {

    // Minimal required interface ----------------------
    Vertex() { }

    void Clear( void * =0 ) {
         point[0] = point[1] = point[2] = 0.0;
    }


    void AddWithWeight(Vertex const & src, double weight) {
        point[0] += weight * src.point[0];
        point[1] += weight * src.point[1];
        point[2] += weight * src.point[2];
    }

    double point[3];
};


// vertex value typedef
struct VertexValue {

	// Basic 'uv' layout channel
	double value;

	// Minimal required interface ----------------------

	void Clear()
	{
		value = 0;
	}


	void AddWithWeight(VertexValue const & src, double weight)
	{
		value += weight * src.value;
	}

};



IndexType GetOSDIndexFromKratosID(IndexType KratosID, std::map<IndexType, IndexType>& rMap)
{
    for (auto it = rMap.begin(); it != rMap.end(); ++it)
        if (it->second == KratosID)
            return it->first;

    // KRATOS_INFO("GetOSDIndexFromKratosID :: KratosID") << KratosID << std::endl;
    // KRATOS_INFO("GetOSDIndexFromKratosID :: rMap loop size") << rMap.size() << std::endl;
    for (auto& it : rMap) {
        // KRATOS_INFO("GetOSDIndexFromKratosID :: osd_index ") << it.first << " kratos id " << it.second << std::endl;

    }

    KRATOS_ERROR << "GetOSDIndexFromKratosID could not find corresponding ID in map!";
}

OpenSubdiv::Far::TopologyRefiner * GenerateOpenSubdivRefiner(
    ModelPart::MeshType& rGrid,
    std::map<IndexType, std::map<IndexType, IndexType>>& rGlobalMapOsdVertToKratosNodeIds,
    std::map<IndexType, std::map<IndexType, IndexType>>& rGlobalMapOsdFaceToKratosCondIds,
    const bool FixFreeEdges = false)
{
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner") << std::endl;
    int num_vertices = rGrid.NumberOfNodes();
    int num_faces = rGrid.NumberOfConditions();

    std::vector<int> num_vertices_per_face(rGrid.NumberOfConditions());
    std::vector<int> vertex_indices_per_face;
    // fill num_vertices_per_face, vertex_indices_per_face
    // loop over conditions
    std::map<IndexType, IndexType> Level0MapOsdVertToKratosNodeIds;
    std::map<IndexType, IndexType> Level0MapOsdFaceToKratosCondIds;

    auto t1_v = high_resolution_clock::now();
    const auto nodes_begin = rGrid.NodesBegin();
#pragma omp parallel for
    for(IndexType vert_index = 0; vert_index < num_vertices; ++vert_index) {
        auto node_it = nodes_begin + vert_index;
#pragma omp critical
        Level0MapOsdVertToKratosNodeIds[vert_index] = node_it->GetId();
    }
    auto t2_v = high_resolution_clock::now();
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: number of vertices ") << num_vertices << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: time needed for loop ") << duration_cast<seconds>(t2_v-t1_v).count() << std::endl;

    auto t1_f = high_resolution_clock::now();
    const auto conds_begin = rGrid.ConditionsBegin();
#pragma omp parallel for ordered
    for(IndexType face_index = 0; face_index < num_faces; ++face_index) {
        auto cond_it = conds_begin + face_index;
        Geometry<NodeType>& r_geom_i = cond_it->GetGeometry();
        int num_vertices = r_geom_i.PointsNumber();

        for (int j = 0; j < num_vertices; j++) {
            IndexType kratos_node_id = r_geom_i.GetPoint(j).GetId();
            IndexType osd_vert_idx = GetOSDIndexFromKratosID(kratos_node_id, Level0MapOsdVertToKratosNodeIds);

            #pragma omp ordered
            vertex_indices_per_face.push_back(osd_vert_idx);
        }

        #pragma omp critical
        {
            num_vertices_per_face[face_index] = num_vertices;
            Level0MapOsdFaceToKratosCondIds[face_index] = cond_it->GetId();
        }
    }
    auto t2_f = high_resolution_clock::now();
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: number of faces ") << num_faces << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: time needed for loop ") << duration_cast<seconds>(t2_f-t1_f).count() << std::endl;
    // KRATOS_ERROR << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner debug 1") << std::endl;

    rGlobalMapOsdVertToKratosNodeIds[0] = Level0MapOsdVertToKratosNodeIds;
    rGlobalMapOsdFaceToKratosCondIds[0] = Level0MapOsdFaceToKratosCondIds;

    OpenSubdiv::Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
    OpenSubdiv::Sdc::Options options;
    if (FixFreeEdges) options.SetVtxBoundaryInterpolation(OpenSubdiv::Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
        else options.SetVtxBoundaryInterpolation(OpenSubdiv::Sdc::Options::VTX_BOUNDARY_NONE);
    // options.SetFVarLinearInterpolation(OpenSubdiv::Sdc::Options::FVAR_LINEAR_NONE);
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner debug 2") << std::endl;

    Descriptor desc;
    desc.numVertices = num_vertices;
    desc.numFaces = num_faces;
    desc.numVertsPerFace = &(num_vertices_per_face[0]);
    desc.vertIndicesPerFace = &(vertex_indices_per_face[0]);
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner num_vertices") << num_vertices << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner num_faces") << num_faces << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner num_vertices_per_face") << num_vertices_per_face << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner vertex_indices_per_face") << vertex_indices_per_face << std::endl;
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner debug 3") << std::endl;

    // Instantiate a FarTopologyRefiner from the descriptor
    auto t1_r = high_resolution_clock::now();
    OpenSubdiv::Far::TopologyRefiner * refiner = OpenSubdiv::Far::TopologyRefinerFactory<Descriptor>::Create(desc, OpenSubdiv::Far::TopologyRefinerFactory<Descriptor>::Options(type, options));
    auto t2_r = high_resolution_clock::now();
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: time needed creating refiner ") << duration_cast<seconds>(t2_r-t1_r).count() << std::endl;

    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner finished successfully") << std::endl;
    return refiner;
};


void GenerateRefinementLevels(
    OpenSubdiv::Far::TopologyRefiner * refiner,
    SizeType n_max_ref = 4
)
{
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateRefinementLevels") << std::endl;

    OpenSubdiv::Far::TopologyRefiner::UniformOptions uniform_ref_options = OpenSubdiv::Far::TopologyRefiner::UniformOptions(n_max_ref);
    // generate all topological data in last refinement, since it is not done by default to save storage space
    uniform_ref_options.fullTopologyInLastLevel = true;
    uniform_ref_options.SetRefinementLevel(n_max_ref);
    // refine uniformly by n_uniformRef
    refiner->RefineUniform(uniform_ref_options);

    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateRefinementLevels end of function") << std::endl;
};

void RefineGeometry(
    ModelPart::MeshType& rGrid,
    ModelPart& rOutput,
    std::map<IndexType, std::map<IndexType, IndexType>>& rGlobalMapOsdVertToKratosNodeIds,
    std::map<IndexType, std::map<IndexType, IndexType>>& rGlobalMapOsdFaceToKratosCondIds,
    SizeType ref_level = 1,
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry") << std::endl;

    // clear output modelpart
    rOutput.Clear();

    if (!refiner) {
        /// generate opensubdiv refiner
        refiner = GenerateOpenSubdivRefiner(rGrid, rGlobalMapOsdVertToKratosNodeIds, rGlobalMapOsdFaceToKratosCondIds, false);
        // KRATOS_INFO("CATMULL_CLARK :: CreateSecondRingNeighbourhoodsOfNodes :: creating new refiner ") << std::endl;
    }
    else {
        // unrefine to enable subsequent refining
        refiner->Unrefine();
    }

    OpenSubdiv::Far::TopologyRefiner::UniformOptions uniform_ref_options = OpenSubdiv::Far::TopologyRefiner::UniformOptions(ref_level);
    // generate all topological data in last refinement, since it is not done by default to save storage space
    uniform_ref_options.fullTopologyInLastLevel = true;
    // uniform_ref_options.SetRefinementLevel(ref_level);
    // refine uniformly by n_uniformRef
    refiner->RefineUniform(uniform_ref_options);

    // Create a buffer to hold the position of the refined verts and
    // local points, then copy the coarse positions at the beginning.
    std::vector<Vertex> verts(refiner->GetNumVerticesTotal());
    //Vertex *verts = (Vertex*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(Vertex) );

    SizeType nverts = rGrid.NumberOfNodes();
    //std::vector<double> verts_coords(nverts*3);
    // IndexType index = 0;
    const auto nodes_begin = rGrid.NodesBegin();
    // for(auto node_it = rGrid.NodesBegin(); node_it != rGrid.NodesEnd(); node_it++) {
    for(int index = 0; index < nverts; ++index) {
        // NodeType& r_node = *node_it;
        auto node_it = nodes_begin + index;
        verts[index].point[0] = node_it->Coordinates()[0];
        verts[index].point[1] = node_it->Coordinates()[1];
        verts[index].point[2] = node_it->Coordinates()[2];
        // index += 1;
    }

    //memcpy(&verts[0], &verts_coords[0], nverts*3*sizeof(double));

    // Interpolate vertex primvar data at each level
    OpenSubdiv::Far::PrimvarRefiner primvarRefiner(*refiner);

    Vertex * src = &verts[0];
    for (SizeType level = 1; level <= ref_level; ++level) {
        Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }
    // return the current vertices and faces
    OpenSubdiv::Far::TopologyLevel const & refLastLevel = refiner->GetLevel(ref_level);

    // create maps for storing relation between Kratos and OSD
    std::map<IndexType, IndexType> LevelMapOsdVertToKratosNodeIds;
    std::map<IndexType, IndexType> LevelMapOsdFaceToKratosCondIds;

    nverts = refLastLevel.GetNumVertices();

    // copy vertex information
    int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
    // verts_coords.resize(nverts*3);
    // memcpy(&verts_coords[0], &verts[firstOfLastVerts], nverts*3*sizeof(double));

    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: before creating nodes with subd data and adding them to rOutput ") << std::endl;
    NodesContainerType::Pointer p_new_nodes_container = make_shared<NodesContainerType>();
    // p_new_nodes_container->reserve(nverts);
    std::vector<double> node_coord(3);
    for (IndexType i = 0; i < nverts; ++i) {
        node_coord[0] = verts[firstOfLastVerts + i].point[0];
        node_coord[1] = verts[firstOfLastVerts + i].point[1];
        node_coord[2] = verts[firstOfLastVerts + i].point[2];
        IndexType osd_vert_idx = i;
        IndexType kratos_node_id = osd_vert_idx + 1;    // to avoid kratos id 0
        NodeTypePointer p_node = make_intrusive<Node>(kratos_node_id, node_coord);
        p_new_nodes_container->push_back(p_node);
        LevelMapOsdVertToKratosNodeIds[osd_vert_idx] = kratos_node_id;
    }

    p_new_nodes_container->Unique();

    rOutput.SetNodes(p_new_nodes_container);
    rOutput.GetCommunicator().LocalMesh().SetNodes(p_new_nodes_container);


    SizeType nfaces = refLastLevel.GetNumFaces();
    std::string condition_name = "SurfaceCondition3D4N";
    for(IndexType osd_face_idx = 0; osd_face_idx < nfaces; osd_face_idx++) {
        IndexType kratos_face_index = osd_face_idx + 1;
        OpenSubdiv::Far::ConstIndexArray face_vertices = refLastLevel.GetFaceVertices(osd_face_idx);
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: debug check 1") << std::endl;
        // PointerVector<Node> points_geometry;
        SizeType number_of_nodes = refLastLevel.GetFaceVertices(osd_face_idx).size();
        NodesArrayType node_array(number_of_nodes);
        std::vector<IndexType> node_indices(number_of_nodes);
        for(SizeType node = 0; node < number_of_nodes; ++node) {
            IndexType osd_node_index = refLastLevel.GetFaceVertices(osd_face_idx)[node];
            IndexType kratos_node_index = osd_node_index + 1;
            node_indices[node] = kratos_node_index;
            // IndexType node_index = refLastLevel.GetFaceVertices(osd_face_idx)[node];
            // NodeTypePointer p_node = p_new_nodes_container->GetContainer().at(node_index);
            // node_array.push_back(p_node);
        }
        Properties::Pointer p_property;
        if (rOutput.HasProperties(0)) {
            p_property = rOutput.pGetProperties(0);
        }
        else {
            p_property = rOutput.CreateNewProperties(0);
        }
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: create condition with id and nodes") << node_indices << std::endl;
        rOutput.CreateNewCondition(condition_name, kratos_face_index, node_indices, p_property);
        LevelMapOsdFaceToKratosCondIds[osd_face_idx] = kratos_face_index;
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: created condition number ") << kratos_face_index << std::endl;
        // auto cond = rOutput.GetCondition(kratos_face_index);
        // auto geom = cond.GetGeometry();
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: new geom nodes ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;
    }
    rOutput.GetCommunicator().LocalMesh().SetConditions(rOutput.pConditions());

    // store maps to corresponding level
    rGlobalMapOsdVertToKratosNodeIds[ref_level] = LevelMapOsdVertToKratosNodeIds;
    rGlobalMapOsdFaceToKratosCondIds[ref_level] = LevelMapOsdFaceToKratosCondIds;

    // for(auto node_it = rOutput.NodesBegin(); node_it != rOutput.NodesEnd(); ++node_it) {
    //     NodeType& r_node = *node_it;
    //     KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: node i ") << r_node.GetId() << std::endl;
    //     KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: coords ") << "[ " << r_node.Coordinates()[0]  << " , " << r_node.Coordinates()[1] << " , " << r_node.Coordinates()[2] << " ] " << std::endl;
    // }

    // for(auto cond_it = rOutput.ConditionsBegin(); cond_it != rOutput.ConditionsEnd(); ++cond_it) {
    //     ConditionType& r_cond = *cond_it;
    //     KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: condition i ") << r_cond.GetId() << std::endl;
    //     auto geom = r_cond.GetGeometry();
    //     KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: nodes ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;
    // }

    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: end of function") << std::endl;

    //free(verts);
};

#endif // OPENSUBDIV_UTILITIES_H define
