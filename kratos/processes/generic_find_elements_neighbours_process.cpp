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


void GenericFindElementalNeighboursProcess::ExecuteInitialize()
{
    KRATOS_TRY

    //finding elemental neighbours of nodes
    FindGlobalNodalElementalNeighboursProcess nodal_neigh_proc(mr_model_part);
    nodal_neigh_proc.Execute();

    //main loop
    block_for_each(mr_model_part.Elements(), [&](Element & rElement) {
        //creating temp face conditions
        Geometry<Node<3> >& geom = rElement.GetGeometry();
        const PointerVector<GeometryType> ElemBoundaries = geom.LocalSpaceDimension()==3 ? geom.GenerateFaces() : geom.GenerateEdges() ;
        const unsigned int nBoundaries = ElemBoundaries.size();
        
        //initializing elem neighbours
        if (rElement.Has(NEIGHBOUR_ELEMENTS)) {
            ElementPointerVector& r_neighbour_elements = rElement.GetValue(NEIGHBOUR_ELEMENTS);
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
            rElemNeighbours(i_boundary) = CheckForNeighbourElems(ElemBoundaries[i_boundary], rElement);
        }
    });

    KRATOS_CATCH("Error finding the elemental neighbours")
}


GlobalPointer<Element> GenericFindElementalNeighboursProcess::CheckForNeighbourElems (const Geometry<Node<3> >& rBoundaryGeom,
                                                 Element & rElement)
{
    //creating a a vector of all the elem pointers
    //therefore the elem of interest will be repeated as many times as nodes in the faces are
    ElementPointerVector PointersOfAllBoundaryNodes;
    const unsigned int nNodes = rBoundaryGeom.size();
    for( unsigned int node_i = 0 ; node_i < nNodes; node_i++){ //starting from node 1
        const ElementPointerVector& rNeighCandidates_i = rBoundaryGeom[node_i].GetValue(NEIGHBOUR_ELEMENTS);
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
        if(PointersOfAllBoundaryNodes[j].Id()!=main_elem_id){ //ignoring main node
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


    //TODO: CheckForNeighbourElems should be done with unordered_map instead
    //unfortunately the operator == is not yet supported by the GlobalPointer, so it does not work
    //when this feature is added, replace the previous code with this one:
    /*
    {
        const int main_elem_id = rElement.Id();
        const ElementPointerVector& rNeighCandidates_0 = rFaceGeom[0].GetValue(NEIGHBOUR_ELEMENTS);
        std::unordered_map< GlobalPointer<Element> , int > CandidatesRepetitions;

        for( unsigned int j = 0 ; j < rNeighCandidates_0.size(); j++){ 
            if(rNeighCandidates_0[j].Id()==main_elem_id){ //removing element from the list
                CandidatesRepetitions.insert(std::make_pair(rNeighCandidates_0(j),1)); //found once
            }
        }
        //now checking in the rest of the geometry nodes if element is repeated.
        for( unsigned int node_i = 1 ; node_i < rFaceGeom.size(); node_i++){ //starting from node 1
            const ElementPointerVector& rNeighCandidates_i = rFaceGeom[node_i].GetValue(NEIGHBOUR_ELEMENTS);
            for( unsigned int candidate_i_k = 0 ; candidate_i_k < rNeighCandidates_i.size(); candidate_i_k++){ //all against all
                auto search = CandidatesRepetitions.find(rNeighCandidates_i(candidate_i_k));
                if (search != CandidatesRepetitions.end()) {
                    search->second++;
                }
            }
        }

        //for(const auto iter = CandidatesRepetitions.begin(); iter != CandidatesRepetitions.end(); ++iter){
        for(const auto iter : CandidatesRepetitions) {
            if(iter.second==rFaceGeom.size()){
                return iter.first;
            }
        }

        //if not found, return pointer to element itself
        return nullptr;
    } */


std::vector<bool> GenericFindElementalNeighboursProcess::HasNeighboursInFaces(const Element& rElement)
{
    std::vector<bool> has_neigh_in_faces;    
    if (rElement.Has(NEIGHBOUR_ELEMENTS)) {
        const ElementPointerVector& r_neighbour_elements = rElement.GetValue(NEIGHBOUR_ELEMENTS);
        has_neigh_in_faces.resize(r_neighbour_elements.size());
        for( unsigned int i = 0 ; i < r_neighbour_elements.size(); i++){ //all against all
            has_neigh_in_faces[i] = (r_neighbour_elements(i).get() != nullptr);
        }
    }
    return has_neigh_in_faces;
}


} /* namespace Kratos.*/
