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


#if !defined(KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED



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
  class FindElementalNeighboursProcess 
	: public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of FindElementalNeighboursProcess
      KRATOS_CLASS_POINTER_DEFINITION(FindElementalNeighboursProcess);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      /// avg_elems ------ expected number of neighbour elements per node., 
      /// avg_nodes ------ expected number of neighbour Nodes 
      /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
      FindElementalNeighboursProcess(ModelPart& model_part, int TDim, unsigned int avg_elems = 10)
		: mr_model_part(model_part)
	{
	mavg_elems = avg_elems;
	mTDim=TDim;
// 	mavg_nodes = avg_nodes;
	}

      /// Destructor.
      virtual ~FindElementalNeighboursProcess()
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
			(in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mavg_elems);
			WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
			rE.erase(rE.begin(),rE.end() );	
		}
		for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
		{
			(ie->GetValue(NEIGHBOUR_ELEMENTS)).reserve(3);
			WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
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
				//KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_CONDITIONS)).size() );
			}
		}	

		//adding the neighbouring conditions to the condition 
		//loop over faces
		if (mTDim==2)
		{
			for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
			{	//face nodes
				Geometry<Node<3> >& geom = (ie)->GetGeometry();
				//vector of the 3 faces around the given face
				(ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
				WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);
				//neighb_face is the vector containing pointers to the three faces around ic
				//neighb_face[0] = neighbour face over edge 1-2 of element ic;
				//neighb_face[1] = neighbour face over edge 2-0 of element ic;
				//neighb_face[2] = neighbour face over edge 0-1 of element ic;
				neighb_elems(0) = CheckForNeighbourElems(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
				neighb_elems(1) = CheckForNeighbourElems(geom[2].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
				neighb_elems(2) = CheckForNeighbourElems(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);

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
		}
		ElementsContainerType& rElems = mr_model_part.Elements();
		for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
		{
			WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
			rE.erase(rE.begin(),rE.end());
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
	  return "FindElementalNeighboursProcess";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "FindElementalNeighboursProcess";
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
	int mTDim;
// 	unsigned int mavg_nodes;
        
        
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

		Element::WeakPointer CheckForNeighbourElems (unsigned int Id_1, unsigned int Id_2, WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
		{	//look for the faces around node Id_1
			for( WeakPointerVector< Element >::iterator i =neighbour_elem.begin(); i != neighbour_elem.end(); i++) 
			{	//look for the nodes of the neighbour faces
				Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
				for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++) 
				{	
					if (neigh_elem_geometry[node_i].Id() == Id_2)
					{
						if(i->Id() != elem->Id())
						{
							return *(i.base());
						}
					}
				}			
			}
			return *(elem.base());
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
      FindElementalNeighboursProcess& operator=(FindElementalNeighboursProcess const& rOther);

      /// Copy constructor.
      //FindElementalNeighboursProcess(FindConditionsNeighboursProcess const& rOther);

        
      ///@}    
        
    }; // Class FindConditionsNeighboursProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    FindElementalNeighboursProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const FindElementalNeighboursProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined 


