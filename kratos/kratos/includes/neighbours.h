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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NEIGHBOURS_H_INCLUDED )
#define  KRATOS_NEIGHBOURS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/weak_pointer_vector.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
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
  template<class TNodeType, class TElementType>
  class Neighbours
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Neighbours
      KRATOS_CLASS_POINTER_DEFINITION(Neighbours);
  
      typedef std::size_t IndexType;

      typedef boost::weak_ptr<TNodeType> NodeWeakPointer;
      
      typedef boost::weak_ptr<TElementType> ElementWeakPointer;
      
      /** An array of pointers to elements. */
      typedef WeakPointerVector<TElementType> NeighbourElementsArrayType;
  
      /** An array of pointers to nodes. */
      typedef WeakPointerVector<TNodeType> NeighbourNodesArrayType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Neighbours(){}

      /// Copy constructor.
      Neighbours(const Neighbours& rOther) 
	: mIndex(rOther.mIndex)
	, mpNeighbourElements(rOther.mpNeighbourElements)
	, mpNeighbourNodes(rOther.mpNeighbourNodes) 
      {
      }
      
      Neighbours(IndexType NewIndex, 
		 typename NeighbourElementsArrayType::Pointer pNewNeighbourElements, 
		 typename NeighbourNodesArrayType::Pointer pNewNeighbourNodes)
	: mIndex(NewIndex)
	, mpNeighbourElements(pNewNeighbourElements)
	, mpNeighbourNodes(pNewNeighbourNodes) 
      {
      }

      Neighbours(IndexType NewIndex) 
	:  mIndex(NewIndex)
      {
      }

      /// Destructor.
      virtual ~Neighbours(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Neighbours& operator=(const Neighbours& rOther)
	{
	  mIndex = rOther.mIndex;
	  mpNeighbourElements = rOther.mpNeighbourElements;
	  mpNeighbourNodes = rOther.mpNeighbourNodes;

	  return *this;
	}
      
      
      ///@}
      ///@name Operations
      ///@{
      
      Neighbours Clone()
      {
	typename NeighbourElementsArrayType::Pointer p_neighbour_elements(new NeighbourElementsArrayType(*mpNeighbourElements)); 
	typename NeighbourNodesArrayType::Pointer p_neighbour_nodes(new NeighbourNodesArrayType(*mpNeighbourNodes));

	return Neighbours(mIndex, p_neighbour_elements, p_neighbour_nodes);
      }
            
      
      ///@}
      ///@name Access
      ///@{ 
      
      IndexType Index()
      {return mIndex;}

      
      const typename NeighbourElementsArrayType::Pointer pNeighbourElements() const
	{
	  return mpNeighbourElements;
	}
      
      typename NeighbourElementsArrayType::Pointer pNeighbourElements()
	{
	  return mpNeighbourElements;
	}
      
      NeighbourElementsArrayType const& NeighbourElements() const
	{
	  return *mpNeighbourElements;
	}
      
      NeighbourElementsArrayType& NeighbourElements()
	{
	  return *mpNeighbourElements;
	}
      
      const typename NeighbourNodesArrayType::Pointer pNeighbourNodes() const
	{
	  return mpNeighbourNodes;
	}
      
      typename NeighbourNodesArrayType::Pointer pNeighbourNodes()
	{
	  return mpNeighbourNodes;
	}
      
      NeighbourNodesArrayType const& NeighbourNodes() const
	{
	  return *mpNeighbourNodes;
	}
      
      NeighbourNodesArrayType& NeighbourNodes()
	{
	  return *mpNeighbourNodes;
	}
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return "Neighbours";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream << "Number of neighbour nodes    : " << mpNeighbourNodes->size() << std::endl;
	  rOStream << "Number of neighbour elements : " << mpNeighbourElements->size() << std::endl;
// 	  rOStream << "Neighbour nodes : " <<;
// 	  for(
// 	  rOStream << "Neighbour elements : " << mpNeighbourElements->size() << std::endl;
	  
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
        
      IndexType mIndex;
        
      typename NeighbourElementsArrayType::Pointer mpNeighbourElements;
        
      typename NeighbourNodesArrayType::Pointer mpNeighbourNodes;
        
        
      ///@} 
      ///@name Private Operators
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
      
      ///@}    
        
    }; // Class Neighbours 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TNodeType, class TElementType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Neighbours<TNodeType, TElementType>& rThis);

  /// output stream function
  template<class TNodeType, class TElementType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Neighbours<TNodeType, TElementType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_H_INCLUDED  defined 


