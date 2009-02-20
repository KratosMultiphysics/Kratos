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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-04-30 08:02:32 $
//   Revision:            $Revision: 1.12 $
//
//


#if !defined(KRATOS_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 
// #include <boost/serialization/base_object.hpp>
// #include <boost/serialization/utility.hpp>
// #include <boost/serialization/list.hpp>



// Project includes
#include "includes/define.h"
// #include "includes/process_info.h"
// #include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"


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
    class Communicator 
    {
    public:
      ///@name  Enum's
      ///@{
      
      
      ///@}
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Communicator
      KRATOS_CLASS_POINTER_DEFINITION(Communicator);
  
      typedef unsigned int IndexType;

      typedef unsigned int SizeType;

      typedef Node<3> NodeType;
      typedef Properties PropertiesType;
      typedef Element ElementType;
      typedef Condition ConditionType;

      typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;

      typedef PointerVector<MeshType> MeshesContainerType;

      /// Nodes container. Which is a vector set of nodes with their Id's as key.
      typedef MeshType::NodesContainerType NodesContainerType;

      /** Iterator over the nodes. This iterator is an indirect
	  iterator over Node::Pointer which turn back a reference to
	  node by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::NodeIterator NodeIterator;

      /** Const iterator over the nodes. This iterator is an indirect
	  iterator over Node::Pointer which turn back a reference to
	  node by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::NodeConstantIterator NodeConstantIterator;

      /** Iterator over the properties. This iterator is an indirect
	  iterator over Properties::Pointer which turn back a reference to
	  properties by * operator and not a pointer for more convenient
	  usage. */

      /// Properties container. Which is a vector set of Properties with their Id's as key.
      typedef MeshType::PropertiesContainerType PropertiesContainerType;

      /** Iterator over the Properties. This iterator is an indirect
	  iterator over Node::Pointer which turn back a reference to
	  node by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::PropertiesIterator PropertiesIterator;

      /** Const iterator over the Properties. This iterator is an indirect
	  iterator over Properties::Pointer which turn back a reference to
	  Properties by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::PropertiesConstantIterator PropertiesConstantIterator;

      /** Iterator over the properties. This iterator is an indirect
	  iterator over Properties::Pointer which turn back a reference to
	  properties by * operator and not a pointer for more convenient
	  usage. */

      /// Element container. A vector set of Elements with their Id's as key.
      typedef MeshType::ElementsContainerType ElementsContainerType;

      /** Iterator over the Elements. This iterator is an indirect
	  iterator over Elements::Pointer which turn back a reference to
	  Element by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::ElementIterator ElementIterator;

      /** Const iterator over the Elements. This iterator is an indirect
	  iterator over Elements::Pointer which turn back a reference to
	  Element by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::ElementConstantIterator ElementConstantIterator;

      /// Condintions container. A vector set of Conditions with their Id's as key.
      typedef MeshType::ConditionsContainerType ConditionsContainerType;

       /** Iterator over the Conditions. This iterator is an indirect
	  iterator over Conditions::Pointer which turn back a reference to
	  Condition by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::ConditionIterator ConditionIterator;

      /** Const iterator over the Conditions. This iterator is an indirect
	  iterator over Conditions::Pointer which turn back a reference to
	  Condition by * operator and not a pointer for more convenient
	  usage. */
      typedef MeshType::ConditionConstantIterator ConditionConstantIterator;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Communicator() 
	{
	  MeshType mesh;
	  mLocalMeshes.push_back(mesh.Clone());
	  mGhostMeshes.push_back(mesh.Clone());
	  mInterfaceMeshes.push_back(mesh.Clone());
	}

       /// Copy constructor.
      Communicator(Communicator const& rOther) 
	: mLocalMeshes(rOther.mLocalMeshes)
	, mGhostMeshes(rOther.mGhostMeshes)
	, mInterfaceMeshes(rOther.mInterfaceMeshes)
      {}
      
     
      /// Destructor.
      virtual ~Communicator()
	{
	}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Communicator& operator=(Communicator const& rOther)
	{
	  mLocalMeshes = rOther.mLocalMeshes;
	  mGhostMeshes = rOther.mGhostMeshes;
	  mInterfaceMeshes = rOther.mInterfaceMeshes;

	  return *this;
	}



      ///@}
      ///@name Access
      ///@{

      
      MeshType::Pointer pLocalMesh(IndexType ThisIndex = 0)
	{
	  return mLocalMeshes(ThisIndex);
	}
      
      MeshType::Pointer pGhostMesh(IndexType ThisIndex = 0)
	{
	  return mGhostMeshes(ThisIndex);
	}
      
      MeshType::Pointer pInterfaceMesh(IndexType ThisIndex = 0)
	{
	  return mInterfaceMeshes(ThisIndex);
	}
      
      const MeshType::Pointer pLocalMesh(IndexType ThisIndex = 0) const
	{
	  return mLocalMeshes(ThisIndex);
	}
      
      const MeshType::Pointer pGhostMesh(IndexType ThisIndex = 0) const
	{
	  return mGhostMeshes(ThisIndex);
	}
      
      const MeshType::Pointer pInterfaceMesh(IndexType ThisIndex = 0) const
	{
	  return mInterfaceMeshes(ThisIndex);
	}
      
      MeshType& LocalMesh(IndexType ThisIndex = 0)
	{
	  return mLocalMeshes[ThisIndex];
	}
      
      MeshType& GhostMesh(IndexType ThisIndex = 0)
	{
	  return mGhostMeshes[ThisIndex];
	}
      
      MeshType& InterfaceMesh(IndexType ThisIndex = 0)
	{
	  return mInterfaceMeshes[ThisIndex];
	}
      
      MeshType const& LocalMesh(IndexType ThisIndex = 0) const
	{
	  return mLocalMeshes[ThisIndex];
	}
      
      MeshType const& GhostMesh(IndexType ThisIndex = 0) const
	{
	  return mGhostMeshes[ThisIndex];
	}
      
      MeshType const& InterfaceMesh(IndexType ThisIndex = 0) const
	{
	  return mInterfaceMeshes[ThisIndex];
	}
      
      MeshesContainerType& LocalMeshes()
	{
	  return mLocalMeshes;
	}
      
      MeshesContainerType& GhostMeshes()
	{
	  return mGhostMeshes;
	}
      
      MeshesContainerType& InterfaceMeshes()
	{
	  return mInterfaceMeshes;
	}
      
      MeshesContainerType const& LocalMeshes() const
	{
	  return mLocalMeshes;
	}
      
      MeshesContainerType const& GhostMeshes() const
	{
	  return mGhostMeshes;
	}
      
      MeshesContainerType const& InterfaceMeshes() const
	{
	  return mInterfaceMeshes;
	}
      
      
      ///@}
      ///@name Operations
      ///@{


      
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
	  return "Communicator";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
		for(IndexType i = 0 ; i <  mLocalMeshes.size() ; i++)
		{
		    rOStream << "    Local Mesh " << i << " : " << std::endl;
		    LocalMesh(i).PrintData(rOStream);
		    rOStream << "    Ghost Mesh " << i << " : " << std::endl;
		    GhostMesh(i).PrintData(rOStream);
		    rOStream << "    Interface Mesh " << i << " : " << std::endl;
		    InterfaceMesh(i).PrintData(rOStream);
		}
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
        

      // To store interfaces local entities
      MeshesContainerType mLocalMeshes;

      // To store interfaces ghost entities
      MeshesContainerType mGhostMeshes;

      // To store interfaces ghost+local entities
      MeshesContainerType mInterfaceMeshes;

      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
//       friend class boost::serialization::access;
      
//       template<class TArchive>
// 	  void serialize(TArchive & ThisArchive, const unsigned int ThisVersion)
// 	  {
// /* 	      ThisArchive & mName & mBufferSize & mCurrentIndex; */
// 	  }

//       void RemoveSolutionStepData(IndexType SolutionStepIndex, MeshType& ThisMesh)
// 	{
// 	  for(NodeIterator i_node = ThisMesh.NodesBegin() ; i_node != ThisMesh.NodesEnd() ; ++i_node)
// 	    i_node->RemoveSolutionStepNodalData(SolutionStepIndex);
// 	}
        
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
      
    }; // Class Communicator 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Communicator& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Communicator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_COMMUNICATOR_H_INCLUDED  defined 


