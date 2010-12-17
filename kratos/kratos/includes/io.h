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
//   Date:                $Date: 2008-10-29 14:31:12 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_IO_H_INCLUDED )
#define  KRATOS_IO_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"


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
  class IO
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of IO
      KRATOS_CLASS_POINTER_DEFINITION(IO);

      typedef Node<3> NodeType;
  
      typedef Mesh<NodeType, Properties, Element, Condition> MeshType;

      typedef MeshType::NodesContainerType NodesContainerType;
  
      typedef MeshType::PropertiesContainerType PropertiesContainerType;
  
      typedef MeshType::ElementsContainerType ElementsContainerType;
  
      typedef MeshType::ConditionsContainerType ConditionsContainerType;
  
      typedef std::vector<std::vector<std::size_t> > ConnectivitiesContainerType;

      typedef std::vector<std::vector<std::size_t> > PartitionIndicesContainerType;
      
      typedef std::vector<std::size_t> PartitionIndicesType;
  
      typedef std::size_t SizeType;
      
      typedef matrix<int> GraphType;

  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      IO(){}

      /// Destructor.
      virtual ~IO(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      virtual bool ReadNode(NodeType& rThisNode)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class.", "")
      }
      
      virtual bool ReadNodes(NodesContainerType& rThisNodes)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void WriteNodes(NodesContainerType const& rThisNodes)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadProperties(Properties& rThisProperties)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadProperties(PropertiesContainerType& rThisProperties)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }
      
      virtual void ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual std::size_t  ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void WriteElements(ElementsContainerType const& rThisElements)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual std::size_t  ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

      virtual void ReadMesh(MeshType & rThisMesh)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void ReadModelPart(ModelPart & rThisModelPart)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
      }

      virtual void DivideInputToPartitions(SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
					    PartitionIndicesType const& NodesPartitions, 
					    PartitionIndicesType const& ElementsPartitions, 
					    PartitionIndicesType const& ConditionsPartitions,
					    PartitionIndicesContainerType const& NodesAllPartitions, 
					    PartitionIndicesContainerType const& ElementsAllPartitions, 
					    PartitionIndicesContainerType const& ConditionsAllPartitions)
      {
        KRATOS_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
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

//       /// Turn back information as a string.
//       virtual std::string Info() const
//  {
//    return "IO";
//  }
      
//       /// Print information about this object.
//       virtual void PrintInfo(std::ostream& rOStream) const
//  {
//    rOStream << "IO";
//  }
      

//       /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const
//  {
//  }
      
            
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
      
      /// Assignment operator.
      IO& operator=(IO const& rOther);

      /// Copy constructor.
      IO(IO const& rOther);

        
      ///@}    
        
    }; // Class IO 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
//                  IO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
//                  const IO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_IO_H_INCLUDED  defined 
