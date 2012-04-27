/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;


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
class FindNodalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindNodalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindNodalNeighboursProcess(ModelPart& model_part, unsigned int avg_elems = 10, unsigned int avg_nodes = 10)
        : mr_model_part(model_part)
    {
        mavg_elems = avg_elems;
        mavg_nodes = avg_nodes;
    }

    /// Destructor.
    virtual ~FindNodalNeighboursProcess()
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

    virtual void Execute()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();

        //first of all the neighbour nodes and elements array are initialized to the guessed size
        //and empties the old entries
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            (in->GetValue(NEIGHBOUR_NODES)).reserve(mavg_nodes);
            WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );

            (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mavg_elems);
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );
        }

        //add the neighbour elements to all the nodes in the mesh
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                //KRATOS_WATCH( pGeom[i] );
                (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
                //KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).size() );
            }
        }

        //adding the neighbouring nodes
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);

            for(unsigned int ie = 0; ie < rE.size(); ie++)
            {
                Element::GeometryType& pGeom = rE[ie].GetGeometry();
                for(unsigned int i = 0; i < pGeom.size(); i++)
                {
                    if(pGeom[i].Id() != in->Id() )
                    {
                        Element::NodeType::WeakPointer temp = pGeom(i);
                        AddUniqueWeakPointer< Node<3> >(in->GetValue(NEIGHBOUR_NODES), temp);
                    }
                }
            }
        }
    }

    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());

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
    virtual std::string Info() const
    {
        return "FindNodalNeighboursProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FindNodalNeighboursProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    unsigned int mavg_elems;
    unsigned int mavg_nodes;


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
    FindNodalNeighboursProcess& operator=(FindNodalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindNodalNeighboursProcess(FindNodalNeighboursProcess const& rOther);


    ///@}

}; // Class FindNodalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined 


