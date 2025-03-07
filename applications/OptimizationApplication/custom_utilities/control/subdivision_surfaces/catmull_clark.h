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

#ifndef CATMULL_CLARK_H
#define CATMULL_CLARK_H

// #pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "geometries/quadrilateral_3d_4.h"
#include "expression/container_expression.h"
// #include "opensubdiv_utilities.h"

// Application includes

namespace Kratos
{

using SizeType = std::size_t;
using IndexType = std::size_t;
using NodeType = Node::NodeType;

///@name Kratos Classes
///@{

class CatmullClarkSDS
{
private:
    // variables
    ModelPart& mrControlPolygon;
    const ModelPart& mrControlledMesh;
    std::map<IndexType, std::vector<IndexType>> mFirstRingNodes;
    std::map<IndexType, std::vector<IndexType>> mFirstRingFaces;
    // OpenSubdiv::Far::TopologyRefiner * mRefiner;

    // functions
    // void CreateFirstRingNeighbourhoodsOfNodes(
    //     ModelPart::MeshType& rInputModelPart, 
    //     std::map<IndexType, std::vector<IndexType>> rOutputMap
    //     );
    // void CreateFirstRingNeighbourhoodsOfElements(
    //     ModelPart::MeshType& rInputModelPart,
    //     std::map<IndexType, std::vector<IndexType>> rOutputMap
    //     );
    SizeType FindNextElement(const GlobalPointersVector<Element>& rElements, Element LastElement, IndexType LastIndex);

public:
    // variables
    

    // functions
    CatmullClarkSDS(ModelPart& rControlPolygon, const ModelPart& mrControlledMesh);

    ~CatmullClarkSDS();


    void CreateMappingMatrix(
        std::vector<double>& rOutputData,
        ModelPart& rControlPolygon,
        const ModelPart& rControlledMesh,
        const bool FixFreeEdges
    );

};

///@}
}

#endif // CATMULL_CLARK_H define
