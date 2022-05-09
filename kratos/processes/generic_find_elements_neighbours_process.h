#if !defined(KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"



namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;
typedef Geometry<Node < 3 > > GeometryType;


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class GenericFindElementalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GenericFindElementalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(GenericFindElementalNeighboursProcess);

    typedef GlobalPointersVector<Element> ElementPointerVector;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenericFindElementalNeighboursProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    ~GenericFindElementalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        KRATOS_TRY

        


        //finding elemental neighbours of nodes
        FindGlobalNodalElementalNeighboursProcess nodal_neigh_proc(mr_model_part);
        nodal_neigh_proc.Execute();
        //prepare NEIGHBOUR_ELEMENTS in each element
        InitializeElementNeighbours();

        //main loop
        block_for_each(mr_model_part.Elements(), [&](Element & rElement) {
            //creating temp face conditions
            Geometry<Node<3> >& geom = rElement.GetGeometry();
            const PointerVector<GeometryType> ElemFaces = geom.GenerateFaces();
            const unsigned int nFaces = ElemFaces.size();
            ElementPointerVector& rElemNeighbours = rElement.GetValue(NEIGHBOUR_ELEMENTS);
            //resizing vector to correct dimension
            rElemNeighbours.resize(nFaces);
            //looping faces of element to find neighb element sharing the same face
            for(unsigned int i_face=0; i_face<nFaces; i_face++){
                rElemNeighbours(i_face) = CheckForNeighbourElems(ElemFaces[i_face], rElement);
            }
        });

        KRATOS_CATCH("Error finding the elemental neighbours")
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindElementalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindElementalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mr_model_part;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeElementNeighbours()
    {
        block_for_each(mr_model_part.Elements(), [&](Element & rElement) {
            const unsigned int nfaces = rElement.GetGeometry().FacesNumber();
            if (rElement.Has(NEIGHBOUR_ELEMENTS)) {
                ElementPointerVector& r_neighbour_elements = rElement.GetValue(NEIGHBOUR_ELEMENTS);
                r_neighbour_elements.reserve(nfaces); 
                r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
            } else {
                ElementPointerVector empty_vector;
                empty_vector.reserve(nfaces); 
                rElement.SetValue(NEIGHBOUR_ELEMENTS, empty_vector);
            }
        });

    }


    GlobalPointer<Element> CheckForNeighbourElems (const Geometry<Node<3> >& rFaceGeom,
                                                 Element & rElement)
    {
        //creating a a vector of all the elem pointers
        //therefore the elem of interest will be repeated as many times as nodes in the faces are
        ElementPointerVector PointersOfAllFaceNodes;
        const unsigned int nNodes = rFaceGeom.size();
        for( unsigned int node_i = 0 ; node_i < nNodes; node_i++){ //starting from node 1
            const ElementPointerVector& rNeighCandidates_i = rFaceGeom[node_i].GetValue(NEIGHBOUR_ELEMENTS);
            for( unsigned int candidate_i_k = 0 ; candidate_i_k < rNeighCandidates_i.size(); candidate_i_k++){ //all against all
                PointersOfAllFaceNodes.push_back(rNeighCandidates_i(candidate_i_k));
            }
        }

        //the element itsef will also be present in the n nodes, so must be avoided
        const unsigned int main_elem_id = rElement.Id();
        //we will take the pointers in the first node as candidates
        const unsigned int nCandidates =  rFaceGeom[0].GetValue(NEIGHBOUR_ELEMENTS).size();
        const unsigned int nCompleteList = PointersOfAllFaceNodes.size();
        for( unsigned int j = 0 ; j < nCandidates; j++){ 
            unsigned int repetitions = 1 ; //the current one must be added
            if(PointersOfAllFaceNodes[j].Id()!=main_elem_id){ //ignoring main node
                for( unsigned int k = nCandidates ; k < nCompleteList; k++){ //starting from pointers belonging to the next node
                    if(&(*PointersOfAllFaceNodes(j))==&(*PointersOfAllFaceNodes(k))){
                        repetitions++;
                    }
                }
                if(repetitions==nNodes){
                    return PointersOfAllFaceNodes(j);
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

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{



    
    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GenericFindElementalNeighboursProcess& operator=(GenericFindElementalNeighboursProcess const& rOther);

    ///@}

}; // Class GenericFindElementalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GenericFindElementalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GenericFindElementalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GENERIC_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined 
