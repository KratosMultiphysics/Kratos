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


#if !defined(KRATOS_MESH_H_INCLUDED )
#define  KRATOS_MESH_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "containers/pointer_vector_map.h"
#include "utilities/indexed_object.h"
#include "geometries/geometry.h"



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
  template<class TNodeType, class TPropertiesType, class TElementType, class TConditionType>
  class Mesh
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Mesh
      KRATOS_CLASS_POINTER_DEFINITION(Mesh);
  
      typedef std::size_t IndexType;

      typedef std::size_t SizeType;

      typedef TNodeType NodeType;

      typedef TPropertiesType PropertiesType;

      typedef Geometry<NodeType> GeometryType;

      typedef TElementType ElementType;

      typedef TConditionType ConditionType;

      /// Nodes container. Which is a vector set of nodes with their Id's as key.
      typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;

      /** Iterator over the nodes. This iterator is an indirect
	  iterator over Node::Pointer which turn back a reference to
	  node by * operator and not a pointer for more convenient
	  usage. */
      typedef typename NodesContainerType::iterator NodeIterator;

      /** Const iterator over the nodes. This iterator is an indirect
	  iterator over Node::Pointer which turn back a reference to
	  node by * operator and not a pointer for more convenient
	  usage. */
      typedef typename NodesContainerType::const_iterator NodeConstantIterator;

      /// Properties container. A vector set of properties with their Id's as key.
      typedef PointerVectorSet<PropertiesType, IndexedObject> PropertiesContainerType;
  
      /** Iterator over the properties. This iterator is an indirect
	  iterator over Properties::Pointer which turn back a reference to
	  properties by * operator and not a pointer for more convenient
	  usage. */
      typedef typename PropertiesContainerType::iterator PropertiesIterator;

      /** Const iterator over the properties. This iterator is an indirect
	  iterator over Properties::Pointer which turn back a reference to
	  properties by * operator and not a pointer for more convenient
	  usage. */
      typedef typename PropertiesContainerType::const_iterator PropertiesConstantIterator;

      /// Geometries container. A vector map of Geometries with given Id's as key.
/*       typedef PointerVectorMap<GeometryType> GeometriesContainerType; */

      /// Element container. A vector set of Elements with their Id's as key.
      typedef PointerVectorSet<ElementType, IndexedObject> ElementsContainerType;

      /** Iterator over the Elements. This iterator is an indirect
	  iterator over Elements::Pointer which turn back a reference to
	  Element by * operator and not a pointer for more convenient
	  usage. */
      typedef typename ElementsContainerType::iterator ElementIterator;

      /** Const iterator over the Elements. This iterator is an indirect
	  iterator over Elements::Pointer which turn back a reference to
	  Element by * operator and not a pointer for more convenient
	  usage. */
      typedef typename ElementsContainerType::const_iterator ElementConstantIterator;

      /// Conditions container. A vector set of Conditions with their Id's as key.
      typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;

       /** Iterator over the Conditions. This iterator is an indirect
	  iterator over Conditions::Pointer which turn back a reference to
	  Condition by * operator and not a pointer for more convenient
	  usage. */
      typedef typename ConditionsContainerType::iterator ConditionIterator;

      /** Const iterator over the Conditions. This iterator is an indirect
	  iterator over Conditions::Pointer which turn back a reference to
	  Condition by * operator and not a pointer for more convenient
	  usage. */
      typedef typename ConditionsContainerType::const_iterator ConditionConstantIterator;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Mesh() : mpNodes(new NodesContainerType())
	, mpProperties(new PropertiesContainerType())
	, mpElements(new ElementsContainerType())
	, mpConditions(new ConditionsContainerType()){}

      /// Copy constructor.
      Mesh(Mesh const& rOther) : mpNodes(rOther.mpNodes)
	, mpProperties(rOther.mpProperties)
	, mpElements(rOther.mpElements)
	, mpConditions(rOther.mpConditions){}

      /// Components constructor.
      Mesh(typename NodesContainerType::Pointer NewNodes, 
	   typename PropertiesContainerType::Pointer NewProperties,  
	   typename ElementsContainerType::Pointer NewElements,  
	   typename ConditionsContainerType::Pointer NewConditions) 
	: mpNodes(NewNodes), mpProperties(NewProperties) , mpElements(NewElements), mpConditions(NewConditions){}


      /// Destructor.
      virtual ~Mesh(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      Mesh Clone()
      {
	typename NodesContainerType::Pointer p_nodes(new NodesContainerType(*mpNodes));
	typename PropertiesContainerType::Pointer p_poroperties(new PropertiesContainerType(*mpProperties));
	typename ElementsContainerType::Pointer p_elements(new ElementsContainerType(*mpElements));
	typename ConditionsContainerType::Pointer p_conditions(new ConditionsContainerType(*mpConditions));

	return Mesh(p_nodes, p_poroperties, p_elements, p_conditions);
      }
      
      ///@}
      ///@name Nodes
      ///@{ 
      
      SizeType NumberOfNodes() const {return mpNodes->size();}
      
      
      /** Inserts a node in the mesh.
      */
      void AddNode(typename NodeType::Pointer pNewNode)
      {mpNodes->insert(mpNodes->begin(), pNewNode);}

      /** Returns the Node::Pointer  corresponding to it's identifier */
      typename NodeType::Pointer pGetNode(IndexType NodeId){return (*mpNodes)(NodeId);}

      /** Returns a reference node corresponding to it's identifier */
      NodeType& GetNode(IndexType NodeId){return (*mpNodes)[NodeId];}

      /** Remove the node with given Id from mesh.
      */
      void RemoveNode(IndexType NodeId)
      {mpNodes->erase(NodeId);}

      /** Remove given node from mesh.
      */
      void RemoveNode(NodeType& ThisNode)
      {mpNodes->erase(ThisNode.Id());}

      /** Remove given node from mesh.
      */
      void RemoveNode(typename NodeType::Pointer pThisNode)
      {mpNodes->erase(pThisNode->Id());}

      NodeIterator NodesBegin() {return mpNodes->begin();}

      NodeConstantIterator NodesBegin() const {return mpNodes->begin();}

      NodeIterator NodesEnd() {return mpNodes->end();}

      NodeConstantIterator NodesEnd() const {return mpNodes->end();}

      NodesContainerType& Nodes(){return *mpNodes;}

	  typename NodesContainerType::Pointer pNodes(){return mpNodes;}

	  void SetNodes(typename NodesContainerType::Pointer pOtherNodes){mpNodes = pOtherNodes;}

      typename NodesContainerType::ContainerType& NodesArray(){return mpNodes->GetContainer();}

      ///@}
      ///@name Properties
      ///@{ 
      
      SizeType NumberOfProperties() const {return mpProperties->size();}
      
      
      /** Inserts a properties in the mesh.
      */
      void AddProperties(typename PropertiesType::Pointer pNewProperties)
      {mpProperties->insert(mpProperties->begin(), pNewProperties);}

      /** Returns the Properties::Pointer  corresponding to it's identifier */
      typename PropertiesType::Pointer pGetProperties(IndexType PropertiesId){return (*mpProperties)(PropertiesId);}

      /** Returns a reference properties corresponding to it's identifier */
      PropertiesType& GetProperties(IndexType PropertiesId){return (*mpProperties)[PropertiesId];}

      /** Remove the properties with given Id from mesh.
      */
      void RemoveProperties(IndexType PropertiesId)
      {mpProperties->erase(PropertiesId);}

      /** Remove given properties from mesh.
      */
      void RemoveProperties(PropertiesType& ThisProperties)
      {mpProperties->erase(ThisProperties.Id());}

      /** Remove given properties from mesh.
      */
      void RemoveProperties(typename PropertiesType::Pointer pThisProperties)
      {mpProperties->erase(pThisProperties->Id());}

      PropertiesIterator PropertiesBegin() {return mpProperties->begin();}

      PropertiesConstantIterator PropertiesBegin() const {return mpProperties->begin();}

      PropertiesIterator PropertiesEnd() {return mpProperties->end();}

      PropertiesConstantIterator PropertiesEnd() const {return mpProperties->end();}

      PropertiesContainerType& Properties(){return *mpProperties;}

	  typename PropertiesContainerType::Pointer pProperties(){return mpProperties;}

	  void SetProperties(typename PropertiesContainerType::Pointer pOtherProperties){mpProperties = pOtherProperties;}

      typename PropertiesContainerType::ContainerType& PropertiesArray(){return mpProperties->GetContainer();}

      ///@}
      ///@name Elements
      ///@{ 
      
      SizeType NumberOfElements() const {return mpElements->size();}
      
      /** Inserts a element in the mesh.
      */
      void AddElement(typename ElementType::Pointer pNewElement)
      {mpElements->insert(mpElements->begin(), pNewElement);}

      /** Returns the Element::Pointer  corresponding to it's identifier */
      typename ElementType::Pointer pGetElement(IndexType ElementId){return (*mpElements)(ElementId);}

      /** Returns a reference element corresponding to it's identifier */
      ElementType& GetElement(IndexType ElementId){return (*mpElements)[ElementId];}

      /** Remove the element with given Id from mesh.
      */
      void RemoveElement(IndexType ElementId)
      {mpElements->erase(ElementId);}

      /** Remove given element from mesh.
      */
      void RemoveElement(ElementType& ThisElement)
      {mpElements->erase(ThisElement.Id());}

      /** Remove given element from mesh.
      */
      void RemoveElement(typename ElementType::Pointer pThisElement)
      {mpElements->erase(pThisElement->Id());}

      ElementIterator ElementsBegin() {return mpElements->begin();}

      ElementConstantIterator ElementsBegin() const {return mpElements->begin();}

      ElementIterator ElementsEnd() {return mpElements->end();}

      ElementConstantIterator ElementsEnd() const {return mpElements->end();}

      ElementsContainerType& Elements(){return *mpElements;}

	  typename ElementsContainerType::Pointer pElements(){return mpElements;}

	  void SetElements(typename ElementsContainerType::Pointer pOtherElements){mpElements = pOtherElements;}

      typename ElementsContainerType::ContainerType& ElementsArray(){return mpElements->GetContainer();}

      
      ///@}
      ///@name Conditions
      ///@{ 
      
      SizeType NumberOfConditions() const {return mpConditions->size();}
      
      /** Inserts a element in the mesh.
      */
      void AddCondition(typename ConditionType::Pointer pNewCondition)
      {mpConditions->insert(mpConditions->begin(), pNewCondition);}

      /** Returns the Condition::Pointer  corresponding to it's identifier */
      typename ConditionType::Pointer pGetCondition(IndexType ConditionId){return (*mpConditions)(ConditionId);}

      /** Returns a reference condition corresponding to it's identifier */
      ConditionType& GetCondition(IndexType ConditionId){return (*mpConditions)[ConditionId];}

      /** Remove the condition with given Id from mesh.
      */
      void RemoveCondition(IndexType ConditionId)
      {mpConditions->erase(ConditionId);}

      /** Remove given condition from mesh.
      */
      void RemoveCondition(ConditionType& ThisCondition)
      {mpConditions->erase(ThisCondition.Id());}

      /** Remove given condition from mesh.
      */
      void RemoveCondition(typename ConditionType::Pointer pThisCondition)
      {mpConditions->erase(pThisCondition->Id());}

      ConditionIterator ConditionsBegin() {return mpConditions->begin();}

      ConditionConstantIterator ConditionsBegin() const {return mpConditions->begin();}

      ConditionIterator ConditionsEnd() {return mpConditions->end();}

      ConditionConstantIterator ConditionsEnd() const {return mpConditions->end();}

      ConditionsContainerType& Conditions(){return *mpConditions;}

	  typename ConditionsContainerType::Pointer pConditions(){return mpConditions;}

	  void SetConditions(typename ConditionsContainerType::Pointer pOtherConditions){mpConditions = pOtherConditions;}

      typename ConditionsContainerType::ContainerType& ConditionsArray(){return mpConditions->GetContainer();}

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
	  return "Mesh";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream << "    Number of Nodes      : " << mpNodes->size() << std::endl;
 	  rOStream << "    Number of Properties : " << mpProperties->size() << std::endl;
	  rOStream << "    Number of Elements   : " << mpElements->size() << std::endl;
	  rOStream << "    Number of Conditions : " << mpConditions->size() << std::endl;
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
        
      typename NodesContainerType::Pointer mpNodes;

      typename PropertiesContainerType::Pointer mpProperties;
        
      typename ElementsContainerType::Pointer mpElements;

      typename ConditionsContainerType::Pointer mpConditions;
        
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
      Mesh& operator=(const Mesh& rOther)
      {
	mpNodes = rOther.mpNodes;
	mpProperties = rOther.mpProperties;
	mpElements = rOther.mpElements;
	mpConditions = rOther.mpConditions;
      }
      
        
      ///@}    
        
    }; // Class Mesh 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TNodeType, class TPropertiesType, class TElementType, class TConditionType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Mesh<TNodeType, TPropertiesType, TElementType, TConditionType>& rThis);

  /// output stream function
  template<class TNodeType, class TPropertiesType, class TElementType, class TConditionType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Mesh<TNodeType, TPropertiesType, TElementType, TConditionType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_MESH_H_INCLUDED  defined 


