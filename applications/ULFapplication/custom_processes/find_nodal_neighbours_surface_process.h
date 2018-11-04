//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Alex Jarauta
//  Co-author  :     Elaf Mahrous


#if !defined(KRATOS_FIND_NODAL_NEIGHBOURS_SURFACE_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_NODAL_NEIGHBOURS_SURFACE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

#include "processes/process.h"
#include "includes/node.h"
// #include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
// typedef  ModelPart::ElementsContainerType ElementsContainerType;
typedef  ModelPart::ConditionsContainerType ConditionsContainerType;


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
class FindNodalNeighboursSurfaceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindNodalNeighboursSurfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalNeighboursSurfaceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_conds ------ expected number of neighbour conditions per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindNodalNeighboursSurfaceProcess(ModelPart& model_part, const int avg_conds = 8, const int avg_nodes = 8)
        : mr_model_part(model_part)
    {
        mavg_conds = avg_conds;
        mavg_nodes = avg_nodes;
    }

    /// Destructor.
    ~FindNodalNeighboursSurfaceProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
	ConditionsContainerType& rConds = mr_model_part.Conditions();

        //first of all the neighbour nodes and conditions array are initialized to the guessed size
        //and empties the old entries
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
        {
            (in->GetValue(NEIGHBOUR_NODES)).reserve(mavg_nodes);
            WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );
	    
            (in->GetValue(NEIGHBOUR_CONDITIONS)).reserve(mavg_conds);
            WeakPointerVector<Condition >& rC = in->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end() );	    
        }

        //add the neighbour conditions to all the nodes in the mesh
        for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ++ic)
        {
            Condition::GeometryType& pGeom = ic->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                //KRATOS_WATCH( pGeom[i] );
                (pGeom[i].GetValue(NEIGHBOUR_CONDITIONS)).push_back( Condition::WeakPointer( *(ic.base()) ) );
                //KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).size() );
            }
        }

        //adding the neighbouring nodes
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
        {
            WeakPointerVector<Condition >& rC = in->GetValue(NEIGHBOUR_CONDITIONS);

            for(unsigned int ic = 0; ic < rC.size(); ic++)
            {
                Condition::GeometryType& pGeom = rC[ic].GetGeometry();
                for(unsigned int i = 0; i < pGeom.size(); i++)
                {
                    if(pGeom[i].Id() != in->Id() )
                    {
//                         Element::NodeType::WeakPointer temp = pGeom(i);
			//NOT SURE ABOUT THIS!!!!!!!!!
                        Condition::NodeType::WeakPointer temp = pGeom(i);			
                        AddUniqueWeakPointer< Node<3> >(in->GetValue(NEIGHBOUR_NODES), temp);
                    }
                }
            }
        }
    }

    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
        {
            WeakPointerVector<Condition >& rC = in->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end());

            WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );
        }
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
        return "FindNodalNeighboursSurfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalNeighboursSurfaceProcess";
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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mr_model_part;
    int mavg_conds;
    int mavg_nodes;


    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

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
    FindNodalNeighboursSurfaceProcess& operator=(FindNodalNeighboursSurfaceProcess const& rOther);

    /// Copy constructor.
    //FindNodalNeighboursSurfaceProcess(FindNodalNeighboursSurfaceProcess const& rOther);


    ///@}

}; // Class FindNodalNeighboursSurfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalNeighboursSurfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalNeighboursSurfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_NEIGHBOURS_SURFACE_PROCESS_H_INCLUDED  defined 


