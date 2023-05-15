//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Pablo Becker
//

/* Sysytem includes */
#include <functional>

/* Project includes */
#include "processes/generic_find_elements_neighbours_process.h"

namespace Kratos
{


void GenericFindElementalNeighboursProcess::Execute()
{
    KRATOS_TRY

    //finding elemental neighbours of nodes
    FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> nodal_neigh_proc(mrModelPart);
    nodal_neigh_proc.Execute();

    const int current_rank = mrModelPart.GetCommunicator().GetDataCommunicator().Rank();

    //main loop
    block_for_each(mrModelPart.Elements(), [&](Element & rElement) {
        //creating temp face conditions
        const auto& r_geom = rElement.GetGeometry();
        const auto elem_boundaries = r_geom.LocalSpaceDimension()==3 ? r_geom.GenerateFaces() : r_geom.GenerateEdges() ;
        const unsigned int nBoundaries = elem_boundaries.size();

        //initializing elem neighbours
        if (rElement.Has(NEIGHBOUR_ELEMENTS)) {
            auto& r_neighbour_elements = rElement.GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.reserve(nBoundaries);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
        } else {
            ElementPointerVector empty_vector;
            empty_vector.reserve(nBoundaries);
            rElement.SetValue(NEIGHBOUR_ELEMENTS, empty_vector);
        }
        ElementPointerVector& rElemNeighbours = rElement.GetValue(NEIGHBOUR_ELEMENTS);
        //resizing vector to correct dimension
        rElemNeighbours.resize(nBoundaries);

        //looping faces of element to find neighb element sharing the same face
        for(unsigned int i_boundary=0; i_boundary<nBoundaries; i_boundary++){
            rElemNeighbours(i_boundary) = CheckForNeighbourElems(elem_boundaries[i_boundary], rElement, current_rank);
        }
    });

    KRATOS_CATCH("Error finding the elemental neighbours")
}

void GenericFindElementalNeighboursProcess::ExecuteInitialize()
{
    KRATOS_WARNING("GenericFindElementalNeighboursProcess") << "'ExecuteInitialize' call is deprecated. Use 'Execute' instead." << std::endl;
    Execute();
}

GlobalPointer<Element> GenericFindElementalNeighboursProcess::CheckForNeighbourElems (const Geometry<Node >& rBoundaryGeom,
                                                                                      Element & rElement,
                                                                                      const int CurrentRank)
{
    //creating a a vector of all the elem pointers
    //therefore the elem of interest will be repeated as many times as nodes in the faces are
    ElementPointerVector PointersOfAllBoundaryNodes;
    const unsigned int nNodes = rBoundaryGeom.size();
    for( unsigned int node_i = 0 ; node_i < nNodes; node_i++){ //starting from node 1
        const auto& rNeighCandidates_i = rBoundaryGeom[node_i].GetValue(NEIGHBOUR_ELEMENTS);
        for( unsigned int candidate_i_k = 0 ; candidate_i_k < rNeighCandidates_i.size(); candidate_i_k++){ //all against all
            PointersOfAllBoundaryNodes.push_back(rNeighCandidates_i(candidate_i_k));
        }
    }

    //the element itsef will also be present in the n nodes, so must be avoided
    const unsigned int main_elem_id = rElement.Id();
    //we will take the pointers in the first node as candidates
    const unsigned int nCandidates =  rBoundaryGeom[0].GetValue(NEIGHBOUR_ELEMENTS).size();
    const unsigned int nCompleteList = PointersOfAllBoundaryNodes.size();
    for( unsigned int j = 0 ; j < nCandidates; j++){
        unsigned int repetitions = 1 ; //the current one must be added
        if(PointersOfAllBoundaryNodes(j).GetRank() != CurrentRank ||
           PointersOfAllBoundaryNodes[j].Id()!=main_elem_id){ //ignoring main node
            for( unsigned int k = nCandidates ; k < nCompleteList; k++){ //starting from pointers belonging to the next node
                if(  PointersOfAllBoundaryNodes(j).get() == PointersOfAllBoundaryNodes(k).get() ){
                    repetitions++;
                }
            }
            if(repetitions==nNodes){
                return PointersOfAllBoundaryNodes(j);
            }
        }
    }

    //if not found, return pointer to element itself
    return nullptr;
}

std::vector<bool> GenericFindElementalNeighboursProcess::HasNeighboursInFaces(const Element& rElement)
{
    std::vector<bool> has_neigh_in_faces;
    if (rElement.Has(NEIGHBOUR_ELEMENTS)) {
        const auto& r_neighbour_elements = rElement.GetValue(NEIGHBOUR_ELEMENTS);
        has_neigh_in_faces.resize(r_neighbour_elements.size());
        for( unsigned int i = 0 ; i < r_neighbour_elements.size(); i++){ //all against all
            has_neigh_in_faces[i] = (r_neighbour_elements(i).get() != nullptr);
        }
    }
    return has_neigh_in_faces;
}


} /* namespace Kratos.*/
