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

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "geometries/quadrilateral_3d_4.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos
{

using SizeType = std::size_t;
using IndexType = std::size_t;

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

// IndexType GetOppositeIndex(ConstIndexArray IndexArray, const IndexType Index)
// {
//     IndexType opposite_vertex;
//     SizeType max_index = IndexArray.size() - 1;
//     if (Index + 2 > max_index) {
//         return IndexArray[Index + 2 - max_index];
//     }
//     else return IndexArray[Index + 2];
// }

// std::vector<IndexType> ReorderIndicesToStartAtValue(const IndexType Value, ConstIndexArray ValueArray)
// {
//     SizeType starting_index = ValueArray.FindIndex(Value);
//     SizeType max_index = ValueArray.size() - 1;
//     std::vector<IndexType> ReturnArray;

//     for (SizeType i = 0; i <= max_index - starting_index; ++i) {
//         ReturnArray.push_back( ValueArray[starting_index + i] );
//     }
//     for (SizeType i = 0; i < starting_index; ++i) {
//         ReturnArray.push_back( ValueArray[i] );
//     }

//     return ReturnArray;
// }

///@}
}

#endif // CATMULL_CLARK_H define
