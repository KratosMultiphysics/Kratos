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
/// OpenSubdiv
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/patchMap.h>

/// This file contains classes for handling OpenSubdiv's types 
/// and functions for applying OpenSubdivs utilities

using namespace Kratos;

using IndexType = std::size_t;
using SizeType = std::size_t;
using NodeType = Node::NodeType;
using ConditionType = ModelPart::ConditionType;
using GeometryType = Geometry<Node>;
using MeshType = ModelPart::MeshType;
using NodesContainerType = ModelPart::NodesContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;
using NodesArrayType = Condition::NodesArrayType;
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



OpenSubdiv::Far::TopologyRefiner * GenerateOpenSubdivRefiner(
    ModelPart::MeshType& rGrid,
    const bool FixFreeEdges = false)
{
    KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner") << std::endl;
    int num_vertices = rGrid.NumberOfNodes();
    int num_faces = rGrid.NumberOfConditions();

    std::vector<int> num_vertices_per_face(rGrid.NumberOfConditions());
    std::vector<int> vertex_indices_per_face;
    // fill num_vertices_per_face, vertex_indices_per_face
    // loop over conditions
    IndexType index = 0;
    for(auto cond_it = rGrid.ConditionsBegin(); cond_it != rGrid.ConditionsEnd(); cond_it++) {
        Condition& r_condition = *cond_it;
        Geometry<NodeType>& geom_i = r_condition.GetGeometry();
        int num_vertices = r_condition.GetGeometry().PointsNumber();
        num_vertices_per_face[index] = num_vertices;
        for (int j = 0; j < num_vertices; j++) {
            vertex_indices_per_face.push_back(geom_i.GetPoint(j).GetId() - 1);
            // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner :: vertex id ") << geom_i.GetPoint(j).GetId() - 1 << std::endl;
        }
        index += 1;
        // ++index;
    }
    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner debug 1") << std::endl;

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
    OpenSubdiv::Far::TopologyRefiner * refiner = OpenSubdiv::Far::TopologyRefinerFactory<Descriptor>::Create(desc, OpenSubdiv::Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

    // KRATOS_INFO("OPENSUBDIV_UTILS :: GenerateOpenSubdivRefiner finished successfully") << std::endl;
    return refiner;
};


void RefineGeometry(
    ModelPart::MeshType& rGrid, 
    ModelPart& rOutput, 
    bool FixFreeEdges = false, 
    SizeType n_uniformRef = 1, 
    OpenSubdiv::Far::TopologyRefiner * refiner = nullptr
    )
{
    KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry") << std::endl;

    // clear output modelpart 
    rOutput.Clear();

    if (!refiner) {
        /// generate opensubdiv refiner if not passed to function
        refiner = GenerateOpenSubdivRefiner(rGrid, FixFreeEdges);
    }
    else {
        // unrefine to enable subsequent refining 
        refiner->Unrefine();
    }

    OpenSubdiv::Far::TopologyRefiner::UniformOptions uniform_ref_options = OpenSubdiv::Far::TopologyRefiner::UniformOptions(n_uniformRef);
    // generate all topological data in last refinement, since it is not done by default to save storage space
    uniform_ref_options.fullTopologyInLastLevel = true;
    uniform_ref_options.SetRefinementLevel(n_uniformRef);
    // refine uniformly by n_uniformRef
    refiner->RefineUniform(uniform_ref_options);

    // Create a buffer to hold the position of the refined verts and
    // local points, then copy the coarse positions at the beginning.
    Vertex *verts = (Vertex*)std::malloc( refiner->GetNumVerticesTotal() *sizeof(Vertex) );
    
    SizeType nverts = rGrid.NumberOfNodes();
    std::vector<double> verts_coords(nverts*3);
    IndexType index = 0;
    for(auto node_it = rGrid.NodesBegin(); node_it != rGrid.NodesEnd(); node_it++) {
        NodeType& r_node = *node_it;
        verts_coords[3*index + 0] = r_node.Coordinates()[0];
        verts_coords[3*index + 1] = r_node.Coordinates()[1];
        verts_coords[3*index + 2] = r_node.Coordinates()[2];
        index += 1;
    }

    memcpy(&verts[0], &verts_coords[0], nverts*3*sizeof(double));

    // Interpolate vertex primvar data at each level
    OpenSubdiv::Far::PrimvarRefiner primvarRefiner(*refiner);

    Vertex * src = verts;
    for (SizeType level = 1; level <= n_uniformRef; ++level) {
        Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
        primvarRefiner.Interpolate(level, src, dst);
        src = dst;
    }
    // return the current vertices and faces
    OpenSubdiv::Far::TopologyLevel const & refLastLevel = refiner->GetLevel(n_uniformRef);
    
    nverts = refLastLevel.GetNumVertices();
    
    // copy vertex information
    int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
    verts_coords.resize(nverts*3);
    memcpy(&verts_coords[0], &verts[firstOfLastVerts], nverts*3*sizeof(double));

    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: before creating nodes with subd data and adding them to rOutput ") << std::endl;
    NodesContainerType::Pointer p_new_nodes_container = make_shared<NodesContainerType>();
    // p_new_nodes_container->reserve(nverts);
    for (IndexType i = 0; i < nverts; ++i) {
        std::vector<double> node_coord(3);
        node_coord[0] = verts_coords[3*i + 0];
        node_coord[1] = verts_coords[3*i + 1];
        node_coord[2] = verts_coords[3*i + 2];
        NodeTypePointer p_node = make_intrusive<Node>(i+1, node_coord);
        p_new_nodes_container->push_back(p_node);
    }

    rOutput.SetNodes(p_new_nodes_container);
    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: new_nodes_container set as Nodes of rOutput") << std::endl;

    // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: before creating conditions with subd data and adding them to rOutput ") << std::endl;

    SizeType nfaces = refLastLevel.GetNumFaces();
    for(IndexType i = 0; i < nfaces; i++) {
        IndexType osd_face_index = i;
        IndexType kratos_face_index = osd_face_index + 1;
        OpenSubdiv::Far::ConstIndexArray face_vertices = refLastLevel.GetFaceVertices(i);
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: debug check 1") << std::endl;
        // PointerVector<Node> points_geometry;
        SizeType number_of_nodes = refLastLevel.GetFaceVertices(i).size();
        NodesArrayType node_array(number_of_nodes);
        std::vector<IndexType> node_indices(number_of_nodes);
        for(SizeType node = 0; node < number_of_nodes; ++node) {
            IndexType osd_node_index = refLastLevel.GetFaceVertices(i)[node];
            IndexType kratos_node_index = osd_node_index + 1;
            node_indices[node] = kratos_node_index;
            // IndexType node_index = refLastLevel.GetFaceVertices(i)[node];
            // NodeTypePointer p_node = p_new_nodes_container->GetContainer().at(node_index);
            // node_array.push_back(p_node);
        }
        std::string condition_name = "SurfaceCondition3D4N";
        Properties::Pointer p_property = rOutput.pGetProperties(0);
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: create condition with id and nodes") << node_indices << std::endl;
        rOutput.CreateNewCondition(condition_name, kratos_face_index, node_indices, p_property);
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: created condition number ") << kratos_face_index << std::endl;
        auto cond = rOutput.GetCondition(kratos_face_index);
        auto geom = cond.GetGeometry();
        // KRATOS_INFO("OPENSUBDIV_UTILS :: RefineGeometry :: new geom nodes ") << "[ " << geom[0].GetId()  << " , " << geom[1].GetId() << " , " << geom[2].GetId() << " , " << geom[3].GetId() << " ] " << std::endl;
    }

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

    free(verts);
};

void GetRegularPatches(
    ModelPart::MeshType& rGrid, 
    std::map<IndexType, std::vector<IndexType>> rFirstRingConditionNeighbourhoods, 
    bool FixFreeEdges = false, 
    SizeType maxPatchLevel = 0
    ) {
    /*
    For now not exactly sure if only regular patches are returned.
    Since OpenSubdiv is used, the output needs to be checked.

    The guess is that, in case the patch is not regular on the coarse 
    grid, the patch is refined adaptively until the patch is regular.
    Thus the returned ids might not be of the original rGrid.
    */
    
    /// generate opensubdiv refiner
    OpenSubdiv::Far::TopologyRefiner * refiner = GenerateOpenSubdivRefiner(rGrid, FixFreeEdges);

    // get refiner
    refiner->RefineUniform(OpenSubdiv::Far::TopologyRefiner::UniformOptions(maxPatchLevel));
    // OpenSubdiv has adaptive refiner to get "normal" subdivision behaviour around extraordinary vertices

    // create patch table
    OpenSubdiv::Far::PatchTableFactory::Options patchOptions(maxPatchLevel);
    patchOptions.endCapType = OpenSubdiv::Far::PatchTableFactory::Options::ENDCAP_BSPLINE_BASIS;

    OpenSubdiv::Far::PatchTable const * patchTable = OpenSubdiv::Far::PatchTableFactory::Create(*refiner, patchOptions);

    // Compute the total number of points we need to evaluate patchtable.
    // we use local points around extraordinary features.
    SizeType nRefinerVertices = refiner->GetNumVerticesTotal();
    SizeType nLocalPoints = patchTable->GetNumLocalPoints();
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchTable->GetNumLocalPoints() ") << patchTable->GetNumLocalPoints() << std::endl;

    // Create a buffer to hold the position of the refined verts and
    // local points, then copy the coarse positions at the beginning.
    SizeType nverts = rGrid.NumberOfNodes();
    std::vector<double> verts_coords(3*nverts);
    int index = 0;
    for(auto node_it = rGrid.NodesBegin(); node_it != rGrid.NodesEnd(); node_it++) {
        NodeType& r_node = *node_it;
        verts_coords[3*index + 0] = r_node.Coordinates()[0];
        verts_coords[3*index + 1] = r_node.Coordinates()[1];
        verts_coords[3*index + 2] = r_node.Coordinates()[2];
        index += 1;
    }
    std::vector<Vertex> verts(nRefinerVertices + nLocalPoints);
    memcpy(&verts[0], &verts_coords[0], nverts*3*sizeof(double));


    // Adaptive refinement to isolate extraordinary vertices, may result in fewer levels than maxIsolation.
    int nRefinedLevels = refiner->GetNumLevels();
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: refiner->GetNumLevels() ") << refiner->GetNumLevels() << std::endl;
    // OpenSubdiv has adaptive refiner to get "normal" subdivision behaviour around extraordinary vertices
    // not sure but i guess that's the idea

    // Evaluate local points from interpolated vertex primvars.
    patchTable->ComputeLocalPointValues(&verts[0], &verts[nRefinerVertices]);

    int numArrays  = (int) patchTable->GetNumPatchArrays();
    int numPatches = (int) patchTable->GetNumPatchesTotal();
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchTable->GetNumPatchArrays ") << patchTable->GetNumPatchArrays() << std::endl;
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchTable->GetNumPatchesTotal ") << patchTable->GetNumPatchesTotal() << std::endl;
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchTable.GetPatchParamTable()[0].GetFaceId() ") << patchTable->GetPatchParamTable()[0].GetFaceId() << std::endl;
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchTable.GetPatchParamTable()[1].GetFaceId() ") << patchTable->GetPatchParamTable()[1].GetFaceId() << std::endl;
    patchTable->print();

    // Create a OpenSubdiv::Far::PatchMap to help locating patches in the table
    OpenSubdiv::Far::PatchMap patchmap(*patchTable);
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchMap created ") << std::endl;

    // refiner->GetLevel(0).PrintTopology();

    int num_faces = refiner->GetLevel(maxPatchLevel).GetNumFaces();

    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: rGrid.NumberOfConditions() ") << rGrid.NumberOfConditions() << std::endl;
    // for(IndexType i = 1; i <= rGrid.NumberOfConditions(); i++) {
    // for(auto cond_it = rGrid.ConditionsBegin(); cond_it != rGrid.ConditionsEnd(); ++cond_it) {
    for(int face_id = 0; face_id < num_faces; ++face_id) {
        // ConditionType& r_condition = *cond_it;
        // IndexType r_condition_id = r_condition.GetId();
        // int face_id = r_condition_id - 1;
        IndexType r_condition_id = face_id + 1; 
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: r_condition_id ") << r_condition_id << std::endl;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: face_id ") << face_id << std::endl;
        double u = 0.5;
        double v = 0.5;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: inloop debug 1 ") << std::endl;
        OpenSubdiv::Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch( face_id, u, v );
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: inloop debug 2 ") << std::endl;
        // get the control variables for this patch
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: handle ") << handle << std::endl;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: handle->arrayIndex ") << handle->arrayIndex << std::endl;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: handle->patchIndex ") << handle->patchIndex << std::endl;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: handle->vertIndex ") << handle->vertIndex << std::endl;
        OpenSubdiv::Far::PatchDescriptor patchDesc = patchTable->GetPatchDescriptor(*handle);
        auto patch_type = patchDesc.GetType();
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchDesc.GetType() ") << patch_type << std::endl;
        int num_ctrl_pts = patchDesc.GetNumControlVertices();
        patchDesc.print();
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: patchDesc.GetNumControlVertices() ") << num_ctrl_pts << std::endl;
        OpenSubdiv::Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: inloop debug 3 ") << std::endl;
        
        // get index vector only if vertices around current face/condition are regular
        std::vector<IndexType> index_vector;
        index_vector.assign(&cvs[0], &cvs[0]+cvs.size());
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: inloop debug 4 ") << std::endl;
        for(IndexType j = 0; j < cvs.size(); j++) {
            // index_vector.push_back(cvs[j]);
            KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: cvs          entries ") << cvs[j] << std::endl;
            KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: index_vector entries ") << index_vector[j] << std::endl;
        }
        rFirstRingConditionNeighbourhoods[r_condition_id] = index_vector;
        KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: inloop debug 5 ") << std::endl;
    }
    KRATOS_INFO("OPENSUBDIV_UTILS :: GetRegularPatches :: end of function ") << std::endl;

};

#endif // OPENSUBDIV_UTILITIES_H define
