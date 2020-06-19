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

#if !defined(KRATOS_MODEL_PART_H_INCLUDED )
#define  KRATOS_MODEL_PART_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/process_info.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "containers/geometry_container.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/communicator.h"
#include "includes/table.h"
#include "containers/pointer_vector_map.h"
#include "containers/pointer_hash_map_set.h"
#include "input_output/logger.h"
#include "includes/kratos_flags.h"
#include "includes/master_slave_constraint.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/variable_data.h"

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

//forward declaring Model to be avoid cross references
class Model;

/// ModelPart class.

/** Detail class definition.
 */
class KRATOS_API(KRATOS_CORE) ModelPart : public DataValueContainer, public Flags
{
    class GetModelPartName : public std::unary_function<const ModelPart* const, std::string>
    {
    public:
        std::string const& operator()(const ModelPart& rModelPart) const
        {
            return rModelPart.Name();
        }
    };
public:
    ///@name  Enum's
    ///@{

    enum OwnershipType
    {
        Kratos_All,
        Kratos_Local,
        Kratos_Ghost,
        Kratos_Ownership_Size
    };

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPart
    //KRATOS_CLASS_POINTER_DEFINITION(ModelPart); //INTENTIONALLY REMOVING DEFINITION - DO NOT UNCOMMENT

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Dof<double> DofType;
    typedef std::vector< DofType::Pointer > DofsVectorType;
    typedef Variable<double> DoubleVariableType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;

//     typedef PointerVectorSet<DofType, SetIdentityFunction<DofType> > DofsArrayType;
    typedef PointerVectorSet<DofType,
                SetIdentityFunction<DofType>,
                std::less<SetIdentityFunction<DofType>::result_type>,
                std::equal_to<SetIdentityFunction<DofType>::result_type>,
                DofType* > DofsArrayType;


    typedef Node < 3 > NodeType;
    typedef Geometry<NodeType> GeometryType;
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

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

    /// The container of the tables. A vector map of the tables.
    typedef PointerVectorMap<SizeType, TableType> TablesContainerType;

    /** Iterator over the Tables. This iterator is an indirect
    iterator over Tables::Pointer which turn back a reference to
    Table by * operator and not a pointer for more convenient
    usage. */
    typedef TablesContainerType::iterator TableIterator;

    /** Const iterator over the Tables. This iterator is an indirect
    iterator over Tables::Pointer which turn back a reference to
    Table by * operator and not a pointer for more convenient
    usage. */
    typedef TablesContainerType::const_iterator TableConstantIterator;
    /**
     *
     */
    /// The container of the constraints
    typedef MeshType::MasterSlaveConstraintType MasterSlaveConstraintType;
    typedef MeshType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;

    /** Iterator over the constraints. This iterator is an indirect
    iterator over MasterSlaveConstraint::Pointer which turn back a reference to
    MasterSlaveConstraint by * operator and not a pointer for more convenient
    usage. */
    typedef MeshType::MasterSlaveConstraintIteratorType MasterSlaveConstraintIteratorType;

    /** Const iterator over the constraints. This iterator is an indirect
    iterator over MasterSlaveConstraint::Pointer which turn back a reference to
    Table by * operator and not a pointer for more convenient
    usage. */
    typedef MeshType::MasterSlaveConstraintConstantIteratorType MasterSlaveConstraintConstantIteratorType;

    /// The Geometry Container.
    /**
    * Contains all geometries, which can be adressed by specific identifiers.
    */
    typedef GeometryContainer<GeometryType> GeometryContainerType;

    /// Geometry Iterator
    typedef typename GeometryContainerType::GeometryIterator GeometryIterator;

    /// Const Geometry Iterator
    typedef typename GeometryContainerType::GeometryConstantIterator GeometryConstantIterator;

    /// Geometry Hash Map Container. Stores with hash of Ids to corresponding geometries.
    typedef typename GeometryContainerType::GeometriesMapType GeometriesMapType;

    /// The container of the sub model parts. A hash table is used.
    /**
    */
    typedef PointerHashMapSet<ModelPart, std::hash< std::string >, GetModelPartName, Kratos::shared_ptr<ModelPart> >  SubModelPartsContainerType;

    /// Iterator over the sub model parts of this model part.
    /**	Note that this iterator only iterates over the next level of
    	sub model parts and does not go through the hierarchy of the
    	sub model parts
    */
    typedef SubModelPartsContainerType::iterator SubModelPartIterator;

    /// Constant iterator over the sub model parts of this model part.
    /**	Note that this iterator only iterates over the next level of
    	sub model parts and does not go through the hierarchy of the
    	sub model parts
    */
    typedef SubModelPartsContainerType::const_iterator SubModelPartConstantIterator;

    ///@}
    ///@name Flags
    ///@{

    KRATOS_DEFINE_LOCAL_FLAG(ALL_ENTITIES);
    KRATOS_DEFINE_LOCAL_FLAG(OVERWRITE_ENTITIES);

    ///@}
    ///@name Life Cycle
    ///@{


    /// Destructor.
    ~ModelPart() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ModelPart & operator=(ModelPart const& rOther) = delete;

    ///@}
    ///@name Solution Steps
    ///@{

    IndexType CreateSolutionStep();

    IndexType CloneSolutionStep();

    // commented due to a bug, Pooyan.
    //       IndexType CreateTimeStep()
    // 	{
    // 	  IndexType new_index = CreateSolutionStep();
    // 	  mProcessInfo.SetAsTimeStepInfo();

    // 	  return new_index;
    // 	}

    IndexType CloneTimeStep();

    IndexType CreateTimeStep(double NewTime);

    IndexType CloneTimeStep(double NewTime);

    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex);

    //this function returns the "Owner" Model
    Model& GetModel()
    {
        return mrModel;
    }

    ///ATTENTION: this function does not touch the coordinates of the nodes.
    ///It just resets the database values to the values at the beginning of the time step
    void ReduceTimeStep(ModelPart& rModelPart, double NewTime);

    ///@}
    ///@name Nodes
    ///@{

    SizeType NumberOfNodes(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NumberOfNodes();
    }

    /** Inserts a node in the current mesh.
     */
    void AddNode(NodeType::Pointer pNewNode, IndexType ThisIndex = 0);

    /** Inserts a list of nodes in a submodelpart provided their Id. Does nothing if applied to the top model part
     */
    void AddNodes(std::vector<IndexType> const& NodeIds, IndexType ThisIndex = 0);

    /** Inserts a list of pointers to nodes
     */
    template<class TIteratorType >
    void AddNodes(TIteratorType nodes_begin,  TIteratorType nodes_end, IndexType ThisIndex = 0)
    {
        KRATOS_TRY
        ModelPart::NodesContainerType  aux;
        ModelPart::NodesContainerType  aux_root; //they may not exist in the root
        ModelPart* root_model_part = &this->GetRootModelPart();

        for(TIteratorType it = nodes_begin; it!=nodes_end; it++)
        {
            auto it_found = root_model_part->Nodes().find(it->Id());
            if(it_found == root_model_part->NodesEnd()) //node does not exist in the top model part
            {
                aux_root.push_back( *(it.base()) ); //node does not exist
                aux.push_back( *(it.base()) );
            }
            else //if it does exist verify it is the same node
            {
                if(&(*it_found) != &(*it))//check if the pointee coincides
                    KRATOS_ERROR << "attempting to add a new node with Id :" << it_found->Id() << ", unfortunately a (different) node with the same Id already exists" << std::endl;
                else
                    aux.push_back( *(it.base()) );
            }
        }

        //now add to the root model part
        for(auto it = aux_root.begin(); it!=aux_root.end(); it++)
            root_model_part->Nodes().push_back( *(it.base()) );
        root_model_part->Nodes().Unique();

        //add to all of the leaves

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Nodes().push_back( *(it.base()) );

            current_part->Nodes().Unique();

            current_part = current_part->GetParentModelPart();
        }

        KRATOS_CATCH("")
    }

    /** Inserts a node in the current mesh.
     */
    NodeType::Pointer CreateNewNode(int Id, double x, double y, double z, VariablesList::Pointer pNewVariablesList, IndexType ThisIndex = 0);

    NodeType::Pointer CreateNewNode(IndexType Id, double x, double y, double z, IndexType ThisIndex = 0);

    NodeType::Pointer CreateNewNode(IndexType Id, double x, double y, double z, double* pThisData, IndexType ThisIndex = 0);

    NodeType::Pointer CreateNewNode(IndexType NodeId, NodeType const& rSourceNode, IndexType ThisIndex = 0);

    void AssignNode(NodeType::Pointer pThisNode, IndexType ThisIndex = 0);

    /** Returns the Node::Pointer  corresponding to it's identifier */
    NodeType::Pointer pGetNode(IndexType NodeId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pGetNode(NodeId);
    }

    /** Returns a reference node corresponding to it's identifier */
    NodeType& GetNode(IndexType NodeId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).GetNode(NodeId);
    }

    const NodeType& GetNode(IndexType NodeId, IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).GetNode(NodeId);
    }

    /** Remove the node with given Id from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveNode(IndexType NodeId, IndexType ThisIndex = 0);

    /** Remove given node from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveNode(NodeType& ThisNode, IndexType ThisIndex = 0);

    /** Remove given node from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveNode(NodeType::Pointer pThisNode, IndexType ThisIndex = 0);

    /** Remove the node with given Id from mesh with ThisIndex in parents and children.
    */
    void RemoveNodeFromAllLevels(IndexType NodeId, IndexType ThisIndex = 0);

    /** Remove given node from current mesh with ThisIndex in parents and children.
    */
    void RemoveNodeFromAllLevels(NodeType& ThisNode, IndexType ThisIndex = 0);

    /** Remove given node from current mesh with ThisIndex in parents and children.
    */
    void RemoveNodeFromAllLevels(NodeType::Pointer pThisNode, IndexType ThisIndex = 0);

    /** erases all nodes identified by "IdentifierFlag" by removing the pointer.
     * Pointers are erased from this level downwards
     * nodes will be automatically destructured
     * when no pointer is left to them
     */
    void RemoveNodes(Flags IdentifierFlag = TO_ERASE);

    /** erases all nodes identified by "IdentifierFlag" by removing the pointer.
     * Pointers will be erase from all levels
     * nodes will be automatically destructured
     * when no pointer is left to them
     */
    void RemoveNodesFromAllLevels(Flags IdentifierFlag = TO_ERASE);

    /** this function gives back the "root" model part, that is the model_part that has no father (non-const version)*/
    ModelPart& GetRootModelPart();

    /** this function gives back the "root" model part, that is the model_part that has no father (const version)*/
    const ModelPart& GetRootModelPart() const;

    NodeIterator NodesBegin(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).NodesBegin();
    }

    NodeConstantIterator NodesBegin(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NodesBegin();
    }

    NodeIterator NodesEnd(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).NodesEnd();
    }

    NodeConstantIterator NodesEnd(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NodesEnd();
    }

    NodesContainerType& Nodes(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).Nodes();
    }

    const NodesContainerType& Nodes(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).Nodes();
    }

    NodesContainerType::Pointer pNodes(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pNodes();
    }

    void SetNodes(NodesContainerType::Pointer pOtherNodes, IndexType ThisIndex = 0)
    {
        GetMesh(ThisIndex).SetNodes(pOtherNodes);
    }

    NodesContainerType::ContainerType& NodesArray(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).NodesArray();
    }

    template<class TDataType>
    void AddNodalSolutionStepVariable(Variable<TDataType> const& ThisVariable)
    {
        if (!HasNodalSolutionStepVariable(ThisVariable))
        {
            // This error prevents memory leaks if variables are being added to a non-empty modelpart
            KRATOS_ERROR_IF((this->GetRootModelPart()).Nodes().size() != 0)
                << "Attempting to add the variable \"" << ThisVariable.Name()
                << "\" to the model part with name \"" << this->Name() << "\" which is not empty" << std::endl;

            mpVariablesList->Add(ThisVariable);
        }
    }

    template<class TDataType>
    bool HasNodalSolutionStepVariable(Variable<TDataType> const& ThisVariable) const
    {
        return mpVariablesList->Has(ThisVariable);
    }

    VariablesList& GetNodalSolutionStepVariablesList()
    {
        return *mpVariablesList;
    }

    VariablesList const& GetNodalSolutionStepVariablesList() const
    {
        return *mpVariablesList;
    }

    VariablesList::Pointer pGetNodalSolutionStepVariablesList()
    {
        return mpVariablesList;
    }

    void SetNodalSolutionStepVariablesList();

    void SetNodalSolutionStepVariablesList(VariablesList::Pointer pNewVariablesList)
    {
        mpVariablesList = pNewVariablesList;
    }

    SizeType GetNodalSolutionStepDataSize()
    {
        return mpVariablesList->DataSize();
    }

    SizeType GetNodalSolutionStepTotalDataSize()
    {
        return mpVariablesList->DataSize() * mBufferSize;
    }

    ///@}
    ///@name Tables
    ///@{

    SizeType NumberOfTables() const
    {
        return mTables.size();
    }

    /** Inserts a Table
     */
    void AddTable(IndexType TableId, TableType::Pointer pNewTable);

    /** Returns the Table::Pointer  corresponding to it's identifier */
    TableType::Pointer pGetTable(IndexType TableId)
    {
        return mTables(TableId);
    }

    /** Returns a reference to Table corresponding to the identifier */
    TableType& GetTable(IndexType TableId)
    {
        return mTables[TableId];
    }

    /** Remove the Table with given Id from current mesh in this modelpart and all its subs.
    */
    void RemoveTable(IndexType TableId);

    /** Remove the Table with given Id from current mesh in parents, itself and all children.
    */
    void RemoveTableFromAllLevels(IndexType TableId);


    TableIterator TablesBegin()
    {
        return mTables.begin();
    }

    TableConstantIterator TablesBegin() const
    {
        return mTables.begin();
    }

    TableIterator TablesEnd()
    {
        return mTables.end();
    }

    TableConstantIterator TablesEnd() const
    {
        return mTables.end();
    }

    TablesContainerType& Tables()
    {
        return mTables;
    }

    TablesContainerType::ContainerType& TablesArray()
    {
        return mTables.GetContainer();
    }

    ///@}
    ///@name MasterSlaveConstraints
    ///@{

    SizeType NumberOfMasterSlaveConstraints(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NumberOfMasterSlaveConstraints();
    }

    MasterSlaveConstraintContainerType& MasterSlaveConstraints(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).MasterSlaveConstraints();
    }

    const MasterSlaveConstraintContainerType& MasterSlaveConstraints(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).MasterSlaveConstraints();
    }

    MasterSlaveConstraintConstantIteratorType  MasterSlaveConstraintsBegin(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).MasterSlaveConstraintsBegin();
    }

    MasterSlaveConstraintConstantIteratorType  MasterSlaveConstraintsEnd(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).MasterSlaveConstraintsEnd();
    }

    MasterSlaveConstraintIteratorType  MasterSlaveConstraintsBegin(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).MasterSlaveConstraintsBegin();
    }

    MasterSlaveConstraintIteratorType  MasterSlaveConstraintsEnd(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).MasterSlaveConstraintsEnd();
    }

    /** Inserts a master-slave constraint in the current modelpart.
     */
    void AddMasterSlaveConstraint(MasterSlaveConstraintType::Pointer pNewMasterSlaveConstraint, IndexType ThisIndex = 0);

    /** Inserts a list of master-slave constraints to a submodelpart provided their Id. Does nothing if applied to the top model part
     */
    void AddMasterSlaveConstraints(std::vector<IndexType> const& MasterSlaveConstraintIds, IndexType ThisIndex = 0);

    /** Inserts a list of pointers to Master-Slave constraints
     */
    template<class TIteratorType >
    void AddMasterSlaveConstraints(TIteratorType constraints_begin,  TIteratorType constraints_end, IndexType ThisIndex = 0)
    {
        KRATOS_TRY
        ModelPart::MasterSlaveConstraintContainerType  aux;
        ModelPart::MasterSlaveConstraintContainerType  aux_root;
        ModelPart* root_model_part = &this->GetRootModelPart();

        for(TIteratorType it = constraints_begin; it!=constraints_end; it++)
        {
            auto it_found = root_model_part->MasterSlaveConstraints().find(it->Id());
            if(it_found == root_model_part->MasterSlaveConstraintsEnd()) //node does not exist in the top model part
            {
                aux_root.push_back( *(it.base()) );
                aux.push_back( *(it.base()) );
            }
            else //if it does exist verify it is the same node
            {
                if(&(*it_found) != &(*it))//check if the pointee coincides
                    KRATOS_ERROR << "attempting to add a new master-slave constraint with Id :" << it_found->Id() << ", unfortunately a (different) master-slave constraint with the same Id already exists" << std::endl;
                else
                    aux.push_back( *(it.base()) );
            }
        }

        for(auto it = aux_root.begin(); it!=aux_root.end(); it++)
                root_model_part->MasterSlaveConstraints().push_back( *(it.base()) );
        root_model_part->MasterSlaveConstraints().Unique();

        //add to all of the leaves

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->MasterSlaveConstraints().push_back( *(it.base()) );

            current_part->MasterSlaveConstraints().Unique();

            current_part = current_part->GetParentModelPart();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Creates a new master-slave constraint in the current modelpart.
     * @todo Replace these 3 functions by one that perfectly forwards arguments, then just define these 3 interfaces on the pybind side
     */
    MasterSlaveConstraint::Pointer CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                    IndexType Id,
                                                                                    DofsVectorType& rMasterDofsVector,
                                                                                    DofsVectorType& rSlaveDofsVector,
                                                                                    const MatrixType& RelationMatrix,
                                                                                    const VectorType& ConstantVector,
                                                                                    IndexType ThisIndex = 0);

    MasterSlaveConstraint::Pointer CreateNewMasterSlaveConstraint(const std::string& ConstraintName,
                                                                                    IndexType Id,
                                                                                    NodeType& rMasterNode,
                                                                                    const DoubleVariableType& rMasterVariable,
                                                                                    NodeType& rSlaveNode,
                                                                                    const DoubleVariableType& rSlaveVariable,
                                                                                    const double Weight,
                                                                                    const double Constant,
                                                                                    IndexType ThisIndex = 0);

    /**
     * @brief Remove the master-slave constraint with given Id from mesh with ThisIndex in this modelpart and all its subs.
     */
    void RemoveMasterSlaveConstraint(IndexType MasterSlaveConstraintId, IndexType ThisIndex = 0);

    /**
     * @brief Remove given master-slave constraint from mesh with ThisIndex in this modelpart and all its subs.
     */
    void RemoveMasterSlaveConstraint(MasterSlaveConstraintType& ThisMasterSlaveConstraint, IndexType ThisIndex = 0);

    /**
     * @brief Remove the master-slave constraint with given Id from mesh with ThisIndex in parents, itself and children.
     */
    void RemoveMasterSlaveConstraintFromAllLevels(IndexType MasterSlaveConstraintId, IndexType ThisIndex = 0);

    /**
     * @brief Remove given master-slave constraint from mesh with ThisIndex in parents, itself and children.
     */
    void RemoveMasterSlaveConstraintFromAllLevels(MasterSlaveConstraintType& ThisMasterSlaveConstraint, IndexType ThisIndex = 0);

    /**
     * @brief It erases all constraints identified by "IdentifierFlag" by removing the pointer.
     * @details Pointers are erased from this level downwards nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the constraints to remove
     */
    void RemoveMasterSlaveConstraints(Flags IdentifierFlag = TO_ERASE);

    /**
     * @brief It erases all constraints identified by "IdentifierFlag" by removing the pointer.
     * @details Pointers will be erase from all levels nodes will be automatically destructured when no pointer is left to them
     * @param IdentifierFlag The flag that identifies the constraints to remove
     */
    void RemoveMasterSlaveConstraintsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

    /** Returns the MasterSlaveConstraint::Pointer  corresponding to it's identifier */
    MasterSlaveConstraintType::Pointer pGetMasterSlaveConstraint(IndexType ConstraintId, IndexType ThisIndex = 0);

    /** Returns a reference MasterSlaveConstraint corresponding to it's identifier */
    MasterSlaveConstraintType& GetMasterSlaveConstraint(IndexType MasterSlaveConstraintId, IndexType ThisIndex = 0);
    /** Returns a const reference MasterSlaveConstraint corresponding to it's identifier */
    const MasterSlaveConstraintType& GetMasterSlaveConstraint(IndexType MasterSlaveConstraintId, IndexType ThisIndex = 0) const ;

    ///@}
    ///@name Properties
    ///@{

    /**
     * @brief Returns the number of properties of the mesh
     * @param ThisIndex The index identifying the mesh
     * @return The number of properties of the mesh
     */
    SizeType NumberOfProperties(IndexType ThisIndex = 0) const;

    /**
     * @brief Inserts a properties in the current mesh.
     * @param pNewProperties The new property pointer to be added
     * @param ThisIndex The index identifying the mesh
     */
    void AddProperties(PropertiesType::Pointer pNewProperties, IndexType ThisIndex = 0);

    /**
     * @brief Returns if the Properties corresponding to it's identifier exists
     * @param PropertiesId The id identifying the property
     * @param ThisIndex The index identifying the mesh
     * @return True if the properties exist, false otherwise
     */
    bool HasProperties(IndexType PropertiesId, IndexType MeshIndex = 0) const;

    /**
     * @brief Returns if the Properties corresponding to it's identifier exists in any of the model parts
     * @param PropertiesId The id identifying the property
     * @param ThisIndex The index identifying the mesh
     * @return True if the properties exist, false otherwise
     */
    bool RecursivelyHasProperties(IndexType PropertiesId, IndexType MeshIndex = 0) const;

    /**
     * @brief Creates a new property in the current mesh
     * @details If the property is already existing it will crash
     * @param PropertiesId The Id of the new property
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The new created properties
     */
    PropertiesType::Pointer CreateNewProperties(IndexType PropertiesId, IndexType MeshIndex = 0);

    /**
     * @brief Returns the Properties::Pointer  corresponding to it's identifier
     * @details If the property is not existing it will return a warning
     * @param PropertiesId The Id of the new property
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (pointer)
     */
    PropertiesType::Pointer pGetProperties(IndexType PropertiesId, IndexType MeshIndex = 0);

    /**
     * @brief Returns the Properties::Pointer  corresponding to it's identifier (const version)
     * @details If the property is not existing it will return a warning
     * @param PropertiesId The Id of the new property
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (pointer)
     */
    PropertiesType::Pointer pGetProperties(IndexType PropertiesId, IndexType MeshIndex = 0) const;

    /**
     * @brief Returns the Properties::Pointer  corresponding to it's identifier
     * @details If the property is not existing it will return a warning
     * @param PropertiesId The Id of the new property
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (reference)
     */
    PropertiesType& GetProperties(IndexType PropertiesId, IndexType MeshIndex = 0);

    /**
     * @brief Returns the Properties::Pointer  corresponding to it's identifier (const version)
     * @details If the property is not existing it will return a warning
     * @param PropertiesId The Id of the new property
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (reference)
     */
    PropertiesType& GetProperties(IndexType PropertiesId, IndexType MeshIndex = 0) const;

    /**
     * @brief Returns if the sub Properties corresponding to it's address exists
     * @param rAddress The text that indicates the structure of subproperties to iterate and found the property of interest
     * @param ThisIndex The index identifying the mesh
     * @return True if the properties exist, false otherwise
     */
    bool HasProperties(
        const std::string& rAddress,
        IndexType MeshIndex = 0
        ) const;

    /**
     * @brief Returns the sub Properties::Pointer  corresponding to it's address
     * @details If the property is not existing it will return a warning
     * @param rAddress The text that indicates the structure of subproperties to iterate and found the property of interest
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (pointer)
     */
    PropertiesType::Pointer pGetProperties(
        const std::string& rAddress,
        IndexType MeshIndex = 0
        );

    /**
     * @brief Returns the sub Properties::Pointer  corresponding to it's address (const version)
     * @details If the property is not existing it will return a warning
     * @param rAddress The text that indicates the structure of subproperties to iterate and found the property of interest
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (pointer)
     */
    const PropertiesType::Pointer pGetProperties(
        const std::string& rAddress,
        IndexType MeshIndex = 0
        ) const;

    /**
     * @brief Returns the sub Properties::Pointer  corresponding to it's address
     * @details If the property is not existing it will return a warning
     * @param rAddress The text that indicates the structure of subproperties to iterate and found the property of interest
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (reference)
     */
    PropertiesType& GetProperties(
        const std::string& rAddress,
        IndexType MeshIndex = 0
        );

    /**
     * @brief Returns the sub Properties::Pointer corresponding to it's address (const version)
     * @details If the property is not existing it will return a warning
     * @param rAddress The text that indicates the structure of subproperties to iterate and found the property of interest
     * @param MeshIndex The Id of the mesh (0 by default)
     * @return The desired properties (reference)
     */
    const PropertiesType& GetProperties(
        const std::string& rAddress,
        IndexType MeshIndex = 0
        ) const;

    /** Remove the Properties with given Id from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveProperties(IndexType PropertiesId, IndexType ThisIndex = 0);

    /** Remove given Properties from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveProperties(PropertiesType& ThisProperties, IndexType ThisIndex = 0);

    /** Remove given Properties from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveProperties(PropertiesType::Pointer pThisProperties, IndexType ThisIndex = 0);

    /** Remove the Properties with given Id from mesh with ThisIndex in parents, itself and children.
    */
    void RemovePropertiesFromAllLevels(IndexType PropertiesId, IndexType ThisIndex = 0);

    /** Remove given Properties from mesh with ThisIndex in parents, itself and children.
    */
    void RemovePropertiesFromAllLevels(PropertiesType& ThisProperties, IndexType ThisIndex = 0);

    /** Remove given Properties from mesh with ThisIndex in parents, itself and children.
    */
    void RemovePropertiesFromAllLevels(PropertiesType::Pointer pThisProperties, IndexType ThisIndex = 0);

    PropertiesIterator PropertiesBegin(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).PropertiesBegin();
    }

    PropertiesConstantIterator PropertiesBegin(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).PropertiesBegin();
    }

    PropertiesIterator PropertiesEnd(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).PropertiesEnd();
    }

    PropertiesConstantIterator PropertiesEnd(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).PropertiesEnd();
    }

    /**
     * temporarily renamed Properties() function because the declaration of Properties()
     * here violates the declaration of Properties() in properties.h
     * (janosch, in agreement with pooyan)
     */
    PropertiesContainerType& rProperties(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).Properties();
    }

    PropertiesContainerType::Pointer pProperties(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pProperties();
    }

    void SetProperties(PropertiesContainerType::Pointer pOtherProperties, IndexType ThisIndex = 0)
    {
        GetMesh(ThisIndex).SetProperties(pOtherProperties);
    }

    PropertiesContainerType::ContainerType& PropertiesArray(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).PropertiesArray();
    }

    ///@}
    ///@name Elements
    ///@{

    SizeType NumberOfElements(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NumberOfElements();
    }

    /** Inserts a element in the current mesh.
     */
    void AddElement(ElementType::Pointer pNewElement, IndexType ThisIndex = 0);

    /** Inserts a list of elements to a submodelpart provided their Id. Does nothing if applied to the top model part
     */
    void AddElements(std::vector<IndexType> const& ElementIds, IndexType ThisIndex = 0);

    /** Inserts a list of pointers to nodes
     */
    template<class TIteratorType >
    void AddElements(TIteratorType elements_begin,  TIteratorType elements_end, IndexType ThisIndex = 0)
    {
        KRATOS_TRY
        ModelPart::ElementsContainerType  aux;
        ModelPart::ElementsContainerType  aux_root;
        ModelPart* root_model_part = &this->GetRootModelPart();

        for(TIteratorType it = elements_begin; it!=elements_end; it++)
        {
            auto it_found = root_model_part->Elements().find(it->Id());
            if(it_found == root_model_part->ElementsEnd()) //node does not exist in the top model part
            {
                aux_root.push_back( *(it.base()) );
                aux.push_back( *(it.base()) );
            }
            else //if it does exist verify it is the same node
            {
                if(&(*it_found) != &(*it))//check if the pointee coincides
                    KRATOS_ERROR << "attempting to add a new element with Id :" << it_found->Id() << ", unfortunately a (different) element with the same Id already exists" << std::endl;
                else
                    aux.push_back( *(it.base()) );
            }
        }

        for(auto it = aux_root.begin(); it!=aux_root.end(); it++)
                root_model_part->Elements().push_back( *(it.base()) );
        root_model_part->Elements().Unique();

        //add to all of the leaves

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Elements().push_back( *(it.base()) );

            current_part->Elements().Unique();

            current_part = current_part->GetParentModelPart();
        }

        KRATOS_CATCH("")
    }

    /** Inserts an element in the current mesh.
     */
    ElementType::Pointer CreateNewElement(std::string ElementName, IndexType Id, std::vector<IndexType> ElementNodeIds, PropertiesType::Pointer pProperties, IndexType ThisIndex = 0);

    /** Inserts an element in the current mesh.
     */
    ElementType::Pointer CreateNewElement(std::string ElementName, IndexType Id, Geometry< Node < 3 > >::PointsArrayType pElementNodes, PropertiesType::Pointer pProperties, IndexType ThisIndex = 0);

    /** Returns the Element::Pointer  corresponding to it's identifier */
    ElementType::Pointer pGetElement(IndexType ElementId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pGetElement(ElementId);
    }

    /** Returns a reference element corresponding to it's identifier */
    ElementType& GetElement(IndexType ElementId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).GetElement(ElementId);
    }

    const ElementType& GetElement(IndexType ElementId, IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).GetElement(ElementId);
    }

    /** Remove the element with given Id from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveElement(IndexType ElementId, IndexType ThisIndex = 0);

    /** Remove given element from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveElement(ElementType& ThisElement, IndexType ThisIndex = 0);

    /** Remove given element from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveElement(ElementType::Pointer pThisElement, IndexType ThisIndex = 0);

    /** Remove the element with given Id from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveElementFromAllLevels(IndexType ElementId, IndexType ThisIndex = 0);

    /** Remove given element from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveElementFromAllLevels(ElementType& ThisElement, IndexType ThisIndex = 0);

    /** Remove given element from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveElementFromAllLevels(ElementType::Pointer pThisElement, IndexType ThisIndex = 0);

    /** erases all elements identified by "IdentifierFlag" by removing the pointer.
         * Pointers are erased from this level downwards
         * nodes will be automatically destructured
         * when no pointer is left to them
         */
    void RemoveElements(Flags IdentifierFlag = TO_ERASE);

    /** erases all elements identified by "IdentifierFlag" by removing the pointer.
     * Pointers will be erase from all levels
     * nodes will be automatically destructured
     * when no pointer is left to them
     */
    void RemoveElementsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

    ElementIterator ElementsBegin(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ElementsBegin();
    }

    ElementConstantIterator ElementsBegin(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).ElementsBegin();
    }

    ElementIterator ElementsEnd(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ElementsEnd();
    }

    ElementConstantIterator ElementsEnd(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).ElementsEnd();
    }

    ElementsContainerType& Elements(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).Elements();
    }

    const ElementsContainerType& Elements(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).Elements();
    }

    ElementsContainerType::Pointer pElements(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pElements();
    }

    void SetElements(ElementsContainerType::Pointer pOtherElements, IndexType ThisIndex = 0)
    {
        GetMesh(ThisIndex).SetElements(pOtherElements);
    }

    ElementsContainerType::ContainerType& ElementsArray(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ElementsArray();
    }

    ///@}
    ///@name Conditions
    ///@{

    SizeType NumberOfConditions(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).NumberOfConditions();
    }

    /** Inserts a condition in the current mesh.
     */
    void AddCondition(ConditionType::Pointer pNewCondition, IndexType ThisIndex = 0);

    /** Inserts a list of conditions to a submodelpart provided their Id. Does nothing if applied to the top model part
     */
    void AddConditions(std::vector<IndexType> const& ConditionIds, IndexType ThisIndex = 0);

    /** Inserts a list of pointers to nodes
     */
    template<class TIteratorType >
    void AddConditions(TIteratorType conditions_begin,  TIteratorType conditions_end, IndexType ThisIndex = 0)
    {
        KRATOS_TRY
        ModelPart::ConditionsContainerType  aux;
        ModelPart::ConditionsContainerType  aux_root;
        ModelPart* root_model_part = &this->GetRootModelPart();

        for(TIteratorType it = conditions_begin; it!=conditions_end; it++)
        {
            auto it_found = root_model_part->Conditions().find(it->Id());
            if(it_found == root_model_part->ConditionsEnd()) //node does not exist in the top model part
            {
                aux.push_back( *(it.base()) );
                aux_root.push_back( *(it.base()) );
            }
            else //if it does exist verify it is the same node
            {
                if(&(*it_found) != &(*it))//check if the pointee coincides
                    KRATOS_ERROR << "attempting to add a new Condition with Id :" << it_found->Id() << ", unfortunately a (different) Condition with the same Id already exists" << std::endl;
                else
                    aux.push_back( *(it.base()) );
            }
        }

        //now add to the root model part
        for(auto it = aux_root.begin(); it!=aux_root.end(); it++)
                root_model_part->Conditions().push_back( *(it.base()) );
        root_model_part->Conditions().Unique();

        //add to all of the leaves

        ModelPart* current_part = this;
        while(current_part->IsSubModelPart())
        {
            for(auto it = aux.begin(); it!=aux.end(); it++)
                current_part->Conditions().push_back( *(it.base()) );

            current_part->Conditions().Unique();

            current_part = current_part->GetParentModelPart();
        }

        KRATOS_CATCH("")
    }

    /** Inserts a condition in the current mesh.
     */
    ConditionType::Pointer CreateNewCondition(std::string ConditionName,
            IndexType Id, std::vector<IndexType> ConditionNodeIds,
            PropertiesType::Pointer pProperties, IndexType ThisIndex = 0);

    /** Inserts a condition in the current mesh.
     */
    ConditionType::Pointer CreateNewCondition(std::string ConditionName,
            IndexType Id, Geometry< Node < 3 > >::PointsArrayType pConditionNodes,
            PropertiesType::Pointer pProperties, IndexType ThisIndex = 0);

    /** Returns the Condition::Pointer  corresponding to it's identifier */
    ConditionType::Pointer pGetCondition(IndexType ConditionId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pGetCondition(ConditionId);
    }

    /** Returns a reference condition corresponding to it's identifier */
    ConditionType& GetCondition(IndexType ConditionId, IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).GetCondition(ConditionId);
    }

    const ConditionType& GetCondition(IndexType ConditionId, IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).GetCondition(ConditionId);
    }

    /**  Remove the condition with given Id from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveCondition(IndexType ConditionId, IndexType ThisIndex = 0);

    /** Remove given condition from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveCondition(ConditionType& ThisCondition, IndexType ThisIndex = 0);

    /** Remove given condition from mesh with ThisIndex in this modelpart and all its subs.
    */
    void RemoveCondition(ConditionType::Pointer pThisCondition, IndexType ThisIndex = 0);

    /**  Remove the condition with given Id from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveConditionFromAllLevels(IndexType ConditionId, IndexType ThisIndex = 0);

    /** Remove given condition from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveConditionFromAllLevels(ConditionType& ThisCondition, IndexType ThisIndex = 0);

    /** Remove given condition from mesh with ThisIndex in parents, itself and children.
    */
    void RemoveConditionFromAllLevels(ConditionType::Pointer pThisCondition, IndexType ThisIndex = 0);

    /** erases all elements identified by "IdentifierFlag" by removing the pointer.
    * Pointers are erased from this level downwards
    * nodes will be automatically destructured
    * when no pointer is left to them
    */
    void RemoveConditions(Flags IdentifierFlag = TO_ERASE);

    /** erases all elements identified by "IdentifierFlag" by removing the pointer.
     * Pointers will be erase from all levels
     * nodes will be automatically destructured
     * when no pointer is left to them
     */
    void RemoveConditionsFromAllLevels(Flags IdentifierFlag = TO_ERASE);

    ConditionIterator ConditionsBegin(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ConditionsBegin();
    }

    ConditionConstantIterator ConditionsBegin(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).ConditionsBegin();
    }

    ConditionIterator ConditionsEnd(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ConditionsEnd();
    }

    ConditionConstantIterator ConditionsEnd(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).ConditionsEnd();
    }

    ConditionsContainerType& Conditions(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).Conditions();
    }

    const ConditionsContainerType& Conditions(IndexType ThisIndex = 0) const
    {
        return GetMesh(ThisIndex).Conditions();
    }

    ConditionsContainerType::Pointer pConditions(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).pConditions();
    }

    void SetConditions(ConditionsContainerType::Pointer pOtherConditions, IndexType ThisIndex = 0)
    {
        GetMesh(ThisIndex).SetConditions(pOtherConditions);
    }

    ConditionsContainerType::ContainerType& ConditionsArray(IndexType ThisIndex = 0)
    {
        return GetMesh(ThisIndex).ConditionsArray();
    }

    ///@}
    ///@name Geometry Container
    ///@{

    SizeType NumberOfGeometries() const
    {
        return mGeometries.NumberOfGeometries();
    }


    /// Adds a geometry to the geometry container.
    void AddGeometry(typename GeometryType::Pointer pNewGeometry);


    /// Returns the Geometry::Pointer corresponding to the Id
    typename GeometryType::Pointer pGetGeometry(IndexType GeometryId) {
        return mGeometries.pGetGeometry(GeometryId);
    }

    /// Returns the const Geometry::Pointer corresponding to the Id
    const typename GeometryType::Pointer pGetGeometry(IndexType GeometryId) const {
        return mGeometries.pGetGeometry(GeometryId);
    }

    /// Returns the Geometry::Pointer corresponding to the name
    typename GeometryType::Pointer pGetGeometry(std::string GeometryName) {
        return mGeometries.pGetGeometry(GeometryName);
    }

    /// Returns the Geometry::Pointer corresponding to the name
    const typename GeometryType::Pointer pGetGeometry(std::string GeometryName) const {
        return mGeometries.pGetGeometry(GeometryName);
    }

    /// Returns a reference geometry corresponding to the id
    GeometryType& GetGeometry(IndexType GeometryId) {
        return mGeometries.GetGeometry(GeometryId);
    }

    /// Returns a const reference geometry corresponding to the id
    const GeometryType& GetGeometry(IndexType GeometryId) const {
        return mGeometries.GetGeometry(GeometryId);
    }

    /// Returns a reference geometry corresponding to the name
    GeometryType& GetGeometry(std::string GeometryName) {
        return mGeometries.GetGeometry(GeometryName);
    }

    /// Returns a const reference geometry corresponding to the name
    const GeometryType& GetGeometry(std::string GeometryName) const {
        return mGeometries.GetGeometry(GeometryName);
    }


    /// Checks if has geometry by id.
    bool HasGeometry(IndexType GeometryId) const {
        return mGeometries.HasGeometry(GeometryId);
    }

    /// Checks if has geometry by name.
    bool HasGeometry(std::string GeometryName) const {
        return mGeometries.HasGeometry(GeometryName);
    }


    /// Removes a geometry by id.
    void RemoveGeometry(IndexType GeometryId);

    /// Removes a geometry by name.
    void RemoveGeometry(std::string GeometryName);

    /// Removes a geometry by id from all root and sub model parts.
    void RemoveGeometryFromAllLevels(IndexType GeometryId);

    /// Removes a geometry by name from all root and sub model parts.
    void RemoveGeometryFromAllLevels(std::string GeometryName);


    /// Begin geometry iterator
    GeometryIterator GeometriesBegin() {
        return mGeometries.GeometriesBegin();
    }

    /// Begin geometry const iterator
    GeometryConstantIterator GeometriesBegin() const {
        return mGeometries.GeometriesBegin();
    }

    /// End geometry iterator
    GeometryIterator GeometriesEnd() {
        return mGeometries.GeometriesEnd();
    }

    /// End geometry const iterator
    GeometryConstantIterator GeometriesEnd() const {
        return mGeometries.GeometriesEnd();
    }


    /// Get geometry map containe
    GeometriesMapType& Geometries()
    {
        return mGeometries.Geometries();
    }

    /// Get geometry map containe
    const GeometriesMapType& Geometries() const
    {
        return mGeometries.Geometries();
    }

    ///@}
    ///@name Sub model parts
    ///@{

    SizeType NumberOfSubModelParts() const
    {
        return mSubModelParts.size();
    }

    /** Creates a new sub model part with given name.
    Does nothing if a sub model part with the same name exist.
    */
    ModelPart& CreateSubModelPart(std::string const& NewSubModelPartName);

    /** Returns a reference to the sub_model part with given string name
    	In debug gives an error if does not exist.
    */
    ModelPart& GetSubModelPart(std::string const& SubModelPartName);

    /** Returns a shared pointer to the sub_model part with given string name
    	In debug gives an error if does not exist.
    */
    ModelPart* pGetSubModelPart(std::string const& SubModelPartName);

    /** Remove a sub modelpart with given name.
    */
    void RemoveSubModelPart(std::string const& ThisSubModelPartName);

    /** Remove given sub model part.
    */
    void RemoveSubModelPart(ModelPart& ThisSubModelPart);

    SubModelPartIterator SubModelPartsBegin()
    {
        return mSubModelParts.begin();
    }

    SubModelPartConstantIterator SubModelPartsBegin() const
    {
        return mSubModelParts.begin();
    }

    SubModelPartIterator SubModelPartsEnd()
    {
        return mSubModelParts.end();
    }

    SubModelPartConstantIterator SubModelPartsEnd() const
    {
        return mSubModelParts.end();
    }

    SubModelPartsContainerType& SubModelParts()
    {
        return mSubModelParts;
    }

    const SubModelPartsContainerType& SubModelParts() const
    {
        return mSubModelParts;
    }

    /** Returns a pointer to the Parent ModelPart
     * Returns a pointer to itself if it is not a SubModelPart
    */
    ModelPart* GetParentModelPart() const;

    /** Returns whether this ModelPart has a SubModelPart with a given name
    */
    bool HasSubModelPart(std::string const& ThisSubModelPartName) const;

    ///@}
    ///@name Access
    ///@{

    ProcessInfo& GetProcessInfo()
    {
        return *mpProcessInfo;
    }

    ProcessInfo const& GetProcessInfo() const
    {
        return *mpProcessInfo;
    }

    ProcessInfo::Pointer pGetProcessInfo()
    {
        return mpProcessInfo;
    }

    const ProcessInfo::Pointer pGetProcessInfo() const
    {
        return mpProcessInfo;
    }

    void SetProcessInfo(ProcessInfo::Pointer pNewProcessInfo)
    {
        mpProcessInfo = pNewProcessInfo;
    }

    void SetProcessInfo(ProcessInfo& NewProcessInfo)
    {
        *mpProcessInfo = NewProcessInfo;
    }

    SizeType NumberOfMeshes()
    {
        return mMeshes.size();
    }

    MeshType::Pointer pGetMesh(IndexType ThisIndex = 0)
    {
        return mMeshes(ThisIndex);
    }

    const MeshType::Pointer pGetMesh(IndexType ThisIndex = 0) const
    {
        return mMeshes(ThisIndex);
    }

    MeshType& GetMesh(IndexType ThisIndex = 0)
    {
        return mMeshes[ThisIndex];
    }

    MeshType const& GetMesh(IndexType ThisIndex = 0) const
    {
        return mMeshes[ThisIndex];
    }

    MeshesContainerType& GetMeshes()
    {
        return mMeshes;
    }

    MeshesContainerType const& GetMeshes() const
    {
        return mMeshes;
    }

    std::string& Name()
    {
        return mName;
    }

    std::string const& Name() const
    {
        return mName;
    }

    Communicator& GetCommunicator()
    {
        return *mpCommunicator;
    }

    Communicator const& GetCommunicator() const
    {
        return *mpCommunicator;
    }

    Communicator::Pointer pGetCommunicator()
    {
        return mpCommunicator;
    }

    void SetCommunicator(Communicator::Pointer pNewCommunicator)
    {
        mpCommunicator = pNewCommunicator;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the full name of the model part (including the parents model parts)
     * @details This is evaluated in a recursive way
     * @return The full name of the model part
     */
    std::string FullName() const
    {
        std::string full_name = this->Name();
        if (this->IsSubModelPart()) {
            full_name = this->GetParentModelPart()->FullName() + "." + full_name;
        }
        return full_name;
    }

    /**
     * @brief This method returns the name list of submodelparts
     * @return A vector conrtaining the list of submodelparts contained
     */
    std::vector<std::string> GetSubModelPartNames();

    /**
     * @brief This method sets the suffer size of the model part database
     * @details Must be called on root model part, otherwise error is thrown
     * @param NewBufferSize The new buffer size to be set
     */
    void SetBufferSize(IndexType NewBufferSize);

    /**
     * @brief This method gets the suffer size of the model part database
     * @return mBufferSize The buffer size
     */
    IndexType GetBufferSize() const
    {
        return mBufferSize;
    }

    /// run input validation
    virtual int Check( ProcessInfo& rCurrentProcessInfo ) const;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    bool IsSubModelPart() const
    {
        return (mpParentModelPart != NULL);
    }

    bool IsDistributed() const
    {
        return mpCommunicator->IsDistributed();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    friend class Model;

    /// Default constructor.
    ModelPart(VariablesList::Pointer pVariableList, Model& rOwnerModel);

    /// Constructor with name
    ModelPart(std::string const& NewName,VariablesList::Pointer pVariableList, Model& rOwnerModel);

    /// Constructor with name and bufferSize
    ModelPart(std::string const& NewName, IndexType NewBufferSize,VariablesList::Pointer pVariableList, Model& rOwnerModel);

    /// Copy constructor.
    ModelPart(ModelPart const& rOther) = delete;


    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    std::string mName; /// The name of the model part

    IndexType mBufferSize; /// The buffers size of the database

    ProcessInfo::Pointer mpProcessInfo; /// The process info instance

    TablesContainerType mTables; /// The tables contained on the model part

    MeshesContainerType mMeshes; /// The container of all meshes

    GeometryContainerType mGeometries; /// The container of geometries

    VariablesList::Pointer mpVariablesList; /// The variable list

    Communicator::Pointer mpCommunicator; /// The communicator

    ModelPart* mpParentModelPart = NULL; /// The parent model part of the current model part

    SubModelPartsContainerType mSubModelParts; /// The container of the submodelparts

    Model& mrModel; /// The model which contains this model part

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method trims a string in the different components to access recursively to any subproperty
     * @param rStringName The given name to be trimmed
     * @return The list of indexes
     */
    std::vector<IndexType> TrimComponentName(const std::string& rStringName) const
    {
        std::vector<IndexType> list_indexes;

        std::stringstream ss(rStringName);
        for (std::string index_string; std::getline(ss, index_string, '.'); ) {
            list_indexes.push_back(std::stoi(index_string));
        }

        KRATOS_ERROR_IF(list_indexes.size() == 0) << "Properties:: Empty list of indexes when reading suproperties" << std::endl;

        return list_indexes;
    }

    /**
     * @brief This method sets the suffer size of the submodelparts belonging to the current model part (recursively)
     * @param NewBufferSize The new buffer size to be set
     */
    void SetBufferSizeSubModelParts(IndexType NewBufferSize);


    void SetParentModelPart(ModelPart* pParentModelPart)
    {
        mpParentModelPart = pParentModelPart;
    }

    template <typename TEntitiesContainerType>
    void AddEntities(TEntitiesContainerType const& Source, TEntitiesContainerType& rDestination, Flags Options)
    {
        //if (Options->Is(ALL_ENTITIES))
        //{
        //	if(Options->Is())
        //}
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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

}; // Class ModelPart

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
KRATOS_API(KRATOS_CORE) inline std::istream & operator >>(std::istream& rIStream,
        ModelPart& rThis)
{
    return rIStream;
}
/// output stream function
KRATOS_API(KRATOS_CORE) inline std::ostream & operator <<(std::ostream& rOStream,
        const ModelPart& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


///@}


} // namespace Kratos.

#endif // KRATOS_MODEL_PART_H_INCLUDED  defined
