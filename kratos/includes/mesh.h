//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
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
#include "containers/flags.h"
#include "containers/data_value_container.h"
#include "includes/master_slave_constraint.h"


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

/// Mesh is the second level of abstraction in the data structure which hold Nodes, Elements and Conditions and their Properties.
/**
 * Mesh is the second level of abstraction in the data structure which hold Nodes, Elements and Conditions and their Properties.
 * In other words, Mesh is a complete pack of all type of entities without any additional data associated with them.
 * So a set of Elements and Conditions with their Nodes and Properties can be grouped together as a Mesh and send to
 * procedures like mesh refinement, material optimization, mesh movement or any other procedure which works on entities
 * without needing additional data for their processes.
*/
template<class TNodeType, class TPropertiesType, class TElementType, class TConditionType>
class Mesh : public DataValueContainer, public Flags
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

    typedef MasterSlaveConstraint MasterSlaveConstraintType;

    typedef Mesh<TNodeType, TPropertiesType, TElementType, TConditionType> MeshType;

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

    // Defining the constraint base type

    /// The container of the constraints
    typedef PointerVectorSet<MasterSlaveConstraintType, IndexedObject> MasterSlaveConstraintContainerType;

    /** Iterator over the constraints. This iterator is an indirect
    iterator over MasterSlaveConstraint::Pointer which turn back a reference to
    MasterSlaveConstraint by * operator and not a pointer for more convenient
    usage. */
    typedef typename MasterSlaveConstraintContainerType::iterator MasterSlaveConstraintIteratorType;

    /** Const iterator over the constraints. This iterator is an indirect
    iterator over MasterSlaveConstraint::Pointer which turn back a reference to
    Table by * operator and not a pointer for more convenient
    usage. */
    typedef typename MasterSlaveConstraintContainerType::const_iterator MasterSlaveConstraintConstantIteratorType;



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Mesh() : Flags()
        , mpNodes(new NodesContainerType())
        , mpProperties(new PropertiesContainerType())
        , mpElements(new ElementsContainerType())
        , mpConditions(new ConditionsContainerType())
        , mpMasterSlaveConstraints(new MasterSlaveConstraintContainerType()) {}

    /// Copy constructor.
    Mesh(Mesh const& rOther) : Flags(rOther)
        , mpNodes(rOther.mpNodes)
        , mpProperties(rOther.mpProperties)
        , mpElements(rOther.mpElements)
        , mpConditions(rOther.mpConditions)
        , mpMasterSlaveConstraints(rOther.mpMasterSlaveConstraints) {}

    /// Components constructor.
    Mesh(typename NodesContainerType::Pointer NewNodes,
         typename PropertiesContainerType::Pointer NewProperties,
         typename ElementsContainerType::Pointer NewElements,
         typename ConditionsContainerType::Pointer NewConditions,
         typename MasterSlaveConstraintContainerType::Pointer NewMasterSlaveConditions)
        : Flags(), mpNodes(NewNodes), mpProperties(NewProperties) , mpElements(NewElements), mpConditions(NewConditions), mpMasterSlaveConstraints(NewMasterSlaveConditions) {}


    /// Destructor.
    ~Mesh() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Mesh Clone()
    {
        typename NodesContainerType::Pointer p_nodes(new NodesContainerType(*mpNodes));
        typename PropertiesContainerType::Pointer p_properties(new PropertiesContainerType(*mpProperties));
        typename ElementsContainerType::Pointer p_elements(new ElementsContainerType(*mpElements));
        typename ConditionsContainerType::Pointer p_conditions(new ConditionsContainerType(*mpConditions));
        typename MasterSlaveConstraintContainerType::Pointer p_master_slave_constraints(new MasterSlaveConstraintContainerType(*mpMasterSlaveConstraints));

        return Mesh(p_nodes, p_properties, p_elements, p_conditions, p_master_slave_constraints);
    }

    void Clear()
    {
        Flags::Clear();
        DataValueContainer::Clear();
        mpNodes->clear();
        mpProperties->clear();
        mpElements->clear();
        mpConditions->clear();
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimensional space of the mesh geometries
	@return SizeType, working space dimension of this geometry.
    */

    SizeType WorkingSpaceDimension() const
    {
      SizeType dimension = 3;

      // NOTE: possible segmentacion fault if a Element or Condition
      // is created using the base class of geometry, then the mpGeometryData
      // of the geometry is a null pointer and has not any mWorkingSpaceDimension
      if(NumberOfElements()!=0)
	dimension = (mpElements->begin())->GetGeometry().WorkingSpaceDimension();
      else if(NumberOfConditions()!=0)
	dimension = (mpConditions->begin())->GetGeometry().WorkingSpaceDimension();
      else if(NumberOfNodes()!=0)
	dimension = (mpNodes->begin())->Dimension();

      return dimension;
    }

    ///@}
    ///@name Nodes
    ///@{

    SizeType NumberOfNodes() const
    {
        return mpNodes->size();
    }


    /** Inserts a node in the mesh.
    */
    void AddNode(typename NodeType::Pointer pNewNode)
    {
        mpNodes->insert(mpNodes->begin(), pNewNode);
    }

    /** Returns the Node::Pointer  corresponding to it's identifier */
    typename NodeType::Pointer pGetNode(IndexType NodeId)
    {
        auto i = mpNodes->find(NodeId);
        KRATOS_ERROR_IF(i == mpNodes->end()) << " node index not found: " << NodeId << ".";
        return *i.base();
    }

    /** Returns a reference node corresponding to it's identifier */
    NodeType& GetNode(IndexType NodeId)
    {
        auto i = mpNodes->find(NodeId);
        KRATOS_ERROR_IF(i == mpNodes->end()) << " node index not found: " << NodeId << ".";
        return *i;
    }

    const NodeType& GetNode(IndexType NodeId) const
    {
        auto i = mpNodes->find(NodeId);
        KRATOS_ERROR_IF(i == mpNodes->end()) << " node index not found: " << NodeId << ".";
        return *i;
    }

    /** Remove the node with given Id from mesh.
    */
    void RemoveNode(IndexType NodeId)
    {
        mpNodes->erase(NodeId);
    }

    /** Remove given node from mesh.
    */
    void RemoveNode(NodeType& ThisNode)
    {
        mpNodes->erase(ThisNode.Id());
    }

    /** Remove given node from mesh.
    */
    void RemoveNode(typename NodeType::Pointer pThisNode)
    {
        mpNodes->erase(pThisNode->Id());
    }

    NodeIterator NodesBegin()
    {
        return mpNodes->begin();
    }

    NodeConstantIterator NodesBegin() const
    {
        return mpNodes->begin();
    }

    NodeIterator NodesEnd()
    {
        return mpNodes->end();
    }

    NodeConstantIterator NodesEnd() const
    {
        return mpNodes->end();
    }

    NodesContainerType& Nodes()
    {
        return *mpNodes;
    }

    const NodesContainerType& Nodes() const
    {
        return *mpNodes;
    }

    typename NodesContainerType::Pointer pNodes()
    {
        return mpNodes;
    }

    void SetNodes(typename NodesContainerType::Pointer pOtherNodes)
    {
        mpNodes = pOtherNodes;
    }

    typename NodesContainerType::ContainerType& NodesArray()
    {
        return mpNodes->GetContainer();
    }

	bool HasNode(IndexType NodeId) const
	{
		return (mpNodes->find(NodeId) != mpNodes->end());
	}

    ///@}
    ///@name Properties
    ///@{

    SizeType NumberOfProperties() const
    {
        return mpProperties->size();
    }


    /** Inserts a properties in the mesh.
    */
    void AddProperties(typename PropertiesType::Pointer pNewProperties)
    {
        mpProperties->insert(mpProperties->begin(), pNewProperties);
    }

    /** Returns the Properties::Pointer  corresponding to it's identifier */
    typename PropertiesType::Pointer pGetProperties(IndexType PropertiesId)
    {
        return (*mpProperties)(PropertiesId);
    }

    /** Returns a reference properties corresponding to it's identifier */
    PropertiesType& GetProperties(IndexType PropertiesId)
    {
        return (*mpProperties)[PropertiesId];
    }

    /** Remove the properties with given Id from mesh.
    */
    void RemoveProperties(IndexType PropertiesId)
    {
        mpProperties->erase(PropertiesId);
    }

    /** Remove given properties from mesh.
    */
    void RemoveProperties(PropertiesType& ThisProperties)
    {
        mpProperties->erase(ThisProperties.Id());
    }

    /** Remove given properties from mesh.
    */
    void RemoveProperties(typename PropertiesType::Pointer pThisProperties)
    {
        mpProperties->erase(pThisProperties->Id());
    }

    PropertiesIterator PropertiesBegin()
    {
        return mpProperties->begin();
    }

    PropertiesConstantIterator PropertiesBegin() const
    {
        return mpProperties->begin();
    }

    PropertiesIterator PropertiesEnd()
    {
        return mpProperties->end();
    }

    PropertiesConstantIterator PropertiesEnd() const
    {
        return mpProperties->end();
    }

    PropertiesContainerType& Properties()
    {
        return *mpProperties;
    }

    const PropertiesContainerType& Properties() const
    {
        return *mpProperties;
    }

    typename PropertiesContainerType::Pointer pProperties()
    {
        return mpProperties;
    }

    void SetProperties(typename PropertiesContainerType::Pointer pOtherProperties)
    {
        mpProperties = pOtherProperties;
    }

    typename PropertiesContainerType::ContainerType& PropertiesArray()
    {
        return mpProperties->GetContainer();
    }

	bool HasProperties(IndexType NodeId) const
	{
		return (mpProperties->find(NodeId) != mpProperties->end());
	}

    ///@}
    ///@name Elements
    ///@{

    SizeType NumberOfElements() const
    {
        return mpElements->size();
    }

    /** Inserts a element in the mesh.
    */
    void AddElement(typename ElementType::Pointer pNewElement)
    {
        mpElements->insert(mpElements->begin(), pNewElement);
    }

    /** Returns the Element::Pointer  corresponding to it's identifier */
    typename ElementType::Pointer pGetElement(IndexType ElementId)
    {
        auto i = mpElements->find(ElementId);
        KRATOS_ERROR_IF(i == mpElements->end()) << " element index not found: " << ElementId << ".";
        return *i.base();
    }

    /** Returns a reference element corresponding to it's identifier */
    ElementType& GetElement(IndexType ElementId)
    {
        auto i = mpElements->find(ElementId);
        KRATOS_ERROR_IF(i == mpElements->end()) << " element index not found: " << ElementId << ".";
        return *i;
    }

    const ElementType& GetElement(IndexType ElementId) const
    {
        auto i = mpElements->find(ElementId);
        KRATOS_ERROR_IF(i == mpElements->end()) << " element index not found: " << ElementId << ".";
        return *i;
    }

    /** Remove the element with given Id from mesh.
    */
    void RemoveElement(IndexType ElementId)
    {
        mpElements->erase(ElementId);
    }

    /** Remove given element from mesh.
    */
    void RemoveElement(ElementType& ThisElement)
    {
        mpElements->erase(ThisElement.Id());
    }

    /** Remove given element from mesh.
    */
    void RemoveElement(typename ElementType::Pointer pThisElement)
    {
        mpElements->erase(pThisElement->Id());
    }

    ElementIterator ElementsBegin()
    {
        return mpElements->begin();
    }

    ElementConstantIterator ElementsBegin() const
    {
        return mpElements->begin();
    }

    ElementIterator ElementsEnd()
    {
        return mpElements->end();
    }

    ElementConstantIterator ElementsEnd() const
    {
        return mpElements->end();
    }

    ElementsContainerType& Elements()
    {
        return *mpElements;
    }

    const ElementsContainerType& Elements() const
    {
        return *mpElements;
    }

    typename ElementsContainerType::Pointer pElements()
    {
        return mpElements;
    }

    void SetElements(typename ElementsContainerType::Pointer pOtherElements)
    {
        mpElements = pOtherElements;
    }

    typename ElementsContainerType::ContainerType& ElementsArray()
    {
        return mpElements->GetContainer();
    }


	bool HasElement(IndexType ElementId) const
	{
		return (mpElements->find(ElementId) != mpElements->end());
	}

    ///@}
    ///@name Conditions
    ///@{

    SizeType NumberOfConditions() const
    {
        return mpConditions->size();
    }

    /** Inserts a condition in the mesh.
    */
    void AddCondition(typename ConditionType::Pointer pNewCondition)
    {
        mpConditions->insert(mpConditions->begin(), pNewCondition);
    }

    /** Returns the Condition::Pointer  corresponding to it's identifier */
    typename ConditionType::Pointer pGetCondition(IndexType ConditionId)
    {
        auto i = mpConditions->find(ConditionId);
        KRATOS_ERROR_IF(i == mpConditions->end()) << " condition index not found: " << ConditionId << ".";
        return *i.base();
    }

    /** Returns a reference condition corresponding to it's identifier */
    ConditionType& GetCondition(IndexType ConditionId)
    {
        auto i = mpConditions->find(ConditionId);
        KRATOS_ERROR_IF(i == mpConditions->end()) << " condition index not found: " << ConditionId << ".";
        return *i;
    }

    const ConditionType& GetCondition(IndexType ConditionId) const
    {
        auto i = mpConditions->find(ConditionId);
        KRATOS_ERROR_IF(i == mpConditions->end()) << " condition index not found: " << ConditionId << ".";
        return *i;
    }

    /** Remove the condition with given Id from mesh.
    */
    void RemoveCondition(IndexType ConditionId)
    {
        mpConditions->erase(ConditionId);
    }

    /** Remove given condition from mesh.
    */
    void RemoveCondition(ConditionType& ThisCondition)
    {
        mpConditions->erase(ThisCondition.Id());
    }

    /** Remove given condition from mesh.
    */
    void RemoveCondition(typename ConditionType::Pointer pThisCondition)
    {
        mpConditions->erase(pThisCondition->Id());
    }

    ConditionIterator ConditionsBegin()
    {
        return mpConditions->begin();
    }

    ConditionConstantIterator ConditionsBegin() const
    {
        return mpConditions->begin();
    }

    ConditionIterator ConditionsEnd()
    {
        return mpConditions->end();
    }

    ConditionConstantIterator ConditionsEnd() const
    {
        return mpConditions->end();
    }

    ConditionsContainerType& Conditions()
    {
        return *mpConditions;
    }

    const ConditionsContainerType& Conditions() const
    {
        return *mpConditions;
    }

    typename ConditionsContainerType::Pointer pConditions()
    {
        return mpConditions;
    }

    void SetConditions(typename ConditionsContainerType::Pointer pOtherConditions)
    {
        mpConditions = pOtherConditions;
    }

    typename ConditionsContainerType::ContainerType& ConditionsArray()
    {
        return mpConditions->GetContainer();
    }

	bool HasCondition(IndexType ConditionId) const
	{
		return (mpConditions->find(ConditionId) != mpConditions->end());
	}


    ///@}
    ///@name MasterSlaveConstraints
    ///@{



    SizeType NumberOfMasterSlaveConstraints() const
    {
        return mpMasterSlaveConstraints->size();
    }

    /** Inserts a master-slave constraint  in the mesh.
    */
    void AddMasterSlaveConstraint(typename MasterSlaveConstraintType::Pointer pNewMasterSlaveConstraint)
    {
        mpMasterSlaveConstraints->insert(mpMasterSlaveConstraints->begin(), pNewMasterSlaveConstraint);
    }

    /** Returns the MasterSlaveConstraint::Pointer  corresponding to it's identifier */
    typename MasterSlaveConstraintType::Pointer pGetMasterSlaveConstraint(IndexType MasterSlaveConstraintId)
    {
        auto i = mpMasterSlaveConstraints->find(MasterSlaveConstraintId);
        KRATOS_ERROR_IF(i == mpMasterSlaveConstraints->end()) << " master-slave constraint index not found: " << MasterSlaveConstraintId << ".";
        return *i.base();
    }

    /** Returns a reference master-slave constraint corresponding to it's identifier */
    MasterSlaveConstraintType& GetMasterSlaveConstraint(IndexType MasterSlaveConstraintId)
    {
        auto i = mpMasterSlaveConstraints->find(MasterSlaveConstraintId);
        KRATOS_ERROR_IF(i == mpMasterSlaveConstraints->end()) << " master-slave constraint index not found: " << MasterSlaveConstraintId << ".";
        return *i;
    }

    const MasterSlaveConstraintType& GetMasterSlaveConstraint(IndexType MasterSlaveConstraintId) const
    {
        auto i = mpMasterSlaveConstraints->find(MasterSlaveConstraintId);
        KRATOS_ERROR_IF(i == mpMasterSlaveConstraints->end()) << " master-slave constraint index not found: " << MasterSlaveConstraintId << ".";
        return *i;
    }

    /** Remove the master-slave constraint with given Id from mesh.
    */
    void RemoveMasterSlaveConstraint(IndexType MasterSlaveConstraintId)
    {
        mpMasterSlaveConstraints->erase(MasterSlaveConstraintId);
    }

    /** Remove given master-slave constraint from mesh.
    */
    void RemoveMasterSlaveConstraint(MasterSlaveConstraintType& ThisMasterSlaveConstraint)
    {
        mpMasterSlaveConstraints->erase(ThisMasterSlaveConstraint.Id());
    }

    /** Remove given master-slave constraint from mesh.
    */
    void RemoveMasterSlaveConstraint(typename MasterSlaveConstraintType::Pointer pThisMasterSlaveConstraint)
    {
        mpMasterSlaveConstraints->erase(pThisMasterSlaveConstraint->Id());
    }

    MasterSlaveConstraintIteratorType MasterSlaveConstraintsBegin()
    {
        return mpMasterSlaveConstraints->begin();
    }

    MasterSlaveConstraintConstantIteratorType MasterSlaveConstraintsBegin() const
    {
        return mpMasterSlaveConstraints->begin();
    }

    MasterSlaveConstraintIteratorType MasterSlaveConstraintsEnd()
    {
        return mpMasterSlaveConstraints->end();
    }

    MasterSlaveConstraintConstantIteratorType MasterSlaveConstraintsEnd() const
    {
        return mpMasterSlaveConstraints->end();
    }

    MasterSlaveConstraintContainerType& MasterSlaveConstraints()
    {
        return *mpMasterSlaveConstraints;
    }

    const MasterSlaveConstraintContainerType& MasterSlaveConstraints() const
    {
        return *mpMasterSlaveConstraints;
    }

    typename MasterSlaveConstraintContainerType::Pointer pMasterSlaveConstraints()
    {
        return mpMasterSlaveConstraints;
    }

    typename MasterSlaveConstraintContainerType::ContainerType& MasterSlaveConstraintsArray()
    {
        return mpMasterSlaveConstraints->GetContainer();
    }


    bool HasMasterSlaveConstraint(IndexType MasterSlaveConstraintId) const
    {
            return (mpMasterSlaveConstraints->find(MasterSlaveConstraintId) != mpMasterSlaveConstraints->end());
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
        return "Mesh";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "    Number of Nodes       : " << mpNodes->size() << std::endl;
        rOStream << "    Number of Properties  : " << mpProperties->size() << std::endl;
        rOStream << "    Number of Elements    : " << mpElements->size() << std::endl;
        rOStream << "    Number of Conditions  : " << mpConditions->size() << std::endl;
        rOStream << "    Number of Constraints : " << mpMasterSlaveConstraints->size() << std::endl;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
    {
        rOStream << PrefixString << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString ) const
    {
        rOStream << PrefixString << "    Number of Nodes       : " << mpNodes->size() << std::endl;
        rOStream << PrefixString << "    Number of Properties  : " << mpProperties->size() << std::endl;
        rOStream << PrefixString << "    Number of Elements    : " << mpElements->size() << std::endl;
        rOStream << PrefixString << "    Number of Conditions  : " << mpConditions->size() << std::endl;
        rOStream << PrefixString << "    Number of Constraints : " << mpMasterSlaveConstraints->size() << std::endl;
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

    typename MasterSlaveConstraintContainerType::Pointer mpMasterSlaveConstraints;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{


    friend class Serializer;


    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Nodes",mpNodes);
        rSerializer.save("Properties",mpProperties);
        rSerializer.save("Elements",mpElements);
        rSerializer.save("Conditions",mpConditions);
        rSerializer.save("Conditions",mpMasterSlaveConstraints);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Nodes",mpNodes);
        rSerializer.load("Properties",mpProperties);
        rSerializer.load("Elements",mpElements);
        rSerializer.load("Conditions",mpConditions);
        rSerializer.load("Conditions",mpMasterSlaveConstraints);
    }


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
        Flags::operator =(rOther);
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


