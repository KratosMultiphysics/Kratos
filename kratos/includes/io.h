//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
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


/**
 * @class IO
 * @ingroup KratosCore
 * @brief IO provides different implementation of input output procedures which can be used to read and write with different formats and characteristics.
 * @details IO provides different implementation of input output procedures which can be used to read and write with different formats and characteristics.
 * An automatic configurable IO module is added to these components providing the complete set of solutions necessary for dealing with multi-disciplinary problems.
 * This IO module uses different component lists to adjust itself when reading and writing new concepts originating from different fields of analysis.
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IO
    KRATOS_CLASS_POINTER_DEFINITION(IO);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( READ );
    KRATOS_DEFINE_LOCAL_FLAG( WRITE );
    KRATOS_DEFINE_LOCAL_FLAG( APPEND );
    KRATOS_DEFINE_LOCAL_FLAG( IGNORE_VARIABLES_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG( SKIP_TIMER );
    KRATOS_DEFINE_LOCAL_FLAG( MESH_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( SCIENTIFIC_PRECISION );

    /// Node type definition
    using NodeType = Node;

    /// Geometry type definition
    using GeometryType = Geometry<NodeType>;

    /// Mesh type definition
    using MeshType = Mesh<NodeType, Properties, Element, Condition>;

    /// Nodes container type within MeshType
    using NodesContainerType = typename MeshType::NodesContainerType;

    /// Properties container type within MeshType
    using PropertiesContainerType = typename MeshType::PropertiesContainerType;

    /// Geometry container type within ModelPart
    using GeometryContainerType = typename ModelPart::GeometryContainerType;

    /// The geometry map type within ModelPart
    using GeometriesMapType = typename ModelPart::GeometriesMapType;

    /// Elements container type within MeshType
    using ElementsContainerType = typename MeshType::ElementsContainerType;

    /// Conditions container type within MeshType
    using ConditionsContainerType = typename MeshType::ConditionsContainerType;

    /// MasterSlaveConstraint container type within MeshType
    using MasterSlaveConstraintContainerType = typename MeshType::MasterSlaveConstraintContainerType;

    /// Connectivities container type
    using ConnectivitiesContainerType = std::vector<std::vector<std::size_t>>;

    /// Partition indices container type
    using PartitionIndicesContainerType = std::vector<std::vector<std::size_t>>;

    /// Partition indices type
    using PartitionIndicesType = std::vector<std::size_t>;

    /// Size type definition
    using SizeType = std::size_t;

    /// Graph type definition
    using GraphType = DenseMatrix<int>;

    // Auxiliary struct containing information about the partitioning of the entities in a ModelPart
    struct PartitioningInfo
    {
        GraphType Graph;
        PartitionIndicesType NodesPartitions;                   // Partition where the Node is local
        PartitionIndicesType ElementsPartitions;                // Partition where the Element is local
        PartitionIndicesType ConditionsPartitions;              // Partition where the Condition is local
        PartitionIndicesType ConstraintsPartitions;             // Partition where the MasterSlaveConstraint is local
        PartitionIndicesType GeometriesPartitions;              // Partition where the Geometry is local
        PartitionIndicesContainerType NodesAllPartitions;       // Partitions, in which the Node is present (local & ghost)
        PartitionIndicesContainerType ElementsAllPartitions;    // Partitions, in which the Element is present (local & ghost)
        PartitionIndicesContainerType ConditionsAllPartitions;  // Partitions, in which the Condition is present (local & ghost)
        PartitionIndicesContainerType ConstraintsAllPartitions; // Partitions, in which the MasterSlaveConstraint is present (local & ghost)
        PartitionIndicesContainerType GeometriesAllPartitions;  // Partitions, in which the Geometry is present (local & ghost)
    };;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IO() = default;

    /// Destructor.
    virtual ~IO() = default;

    /// Copy constructor.
    IO(IO const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IO& operator=(IO const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method reads one node
     * @param rThisNode The node to be read
     */
    virtual bool ReadNode(NodeType& rThisNode)
    {
        KRATOS_ERROR << "Calling base class method (ReadNode). Please check the definition of derived class." << std::endl;
    }

    /**
     * @brief This method reads the nodes from an array of nodes
     * @param rThisNodes The array of nodes to be read
     */
    virtual bool ReadNodes(NodesContainerType& rThisNodes)
    {
        KRATOS_ERROR << "Calling base class method (ReadNodes). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the number of nodes
     * @return The number of nodes
     */
    virtual std::size_t ReadNodesNumber()
    {
        KRATOS_ERROR << "Calling base class method (ReadNodesNumber). Please check the definition of derived class." << std::endl;;
    }

    /**
     * @brief This method writes the nodes from an array of nodes
     * @param rThisNodes The array of nodes to be written
     */
    virtual void WriteNodes(NodesContainerType const& rThisNodes)
    {
        KRATOS_ERROR << "Calling base class method (WriteNodes). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads one Properties
     * @param rThisProperties The Properties to be read
     */
    virtual void ReadProperties(Properties& rThisProperties)
    {
        KRATOS_ERROR << "Calling base class method (ReadProperties). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the Properties from an array of Properties
     * @param rThisProperties The array of Properties to be read
     */
    virtual void ReadProperties(PropertiesContainerType& rThisProperties)
    {
        KRATOS_ERROR << "Calling base class method (ReadProperties). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes one Properties
     * @param rThisProperties The Properties to be written
     */
    virtual void WriteProperties(Properties const& rThisProperties)
    {
        KRATOS_ERROR << "Calling base class method (WriteProperties). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the Properties from an array of Properties
     * @param rThisProperties The array of Properties to be written
     */
    virtual void WriteProperties(PropertiesContainerType const& rThisProperties)
    {
        KRATOS_ERROR << "Calling base class method (WriteProperties). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads one geometry
     * @param rThisNodes The nodes constituting the geometry
     * @param pThisGeometries The pointer to the geometry
     */
    virtual void ReadGeometry(
        NodesContainerType& rThisNodes,
        GeometryType::Pointer& pThisGeometry
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadGeometry). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads an array of geometries
     * @param rThisNodes The nodes constituting the geometry
     * @param rThisGeometry The array of geometries
     */
    virtual void ReadGeometries(
        NodesContainerType& rThisNodes,
        GeometryContainerType& rThisGeometries
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadGeometries). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the geometries connectivities
     * @param rGeometriesConnectivities The geometries connectivities
     * @return The number of geometries
     */
    virtual std::size_t ReadGeometriesConnectivities(ConnectivitiesContainerType& rGeometriesConnectivities)
    {
        KRATOS_ERROR << "Calling base class method (ReadGeometriesConnectivities). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes an array of geometries
     * @param rThisGeometries The array of geometries to be written
     */
    virtual void WriteGeometries(GeometryContainerType const& rThisGeometries)
    {
        KRATOS_ERROR << "Calling base class method (WriteGeometries). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads one element
     * @param rThisNodes The nodes constituying the element
     * @param rThisProperties The Properties of the element
     * @param pThisElements The pointer to the element
     */
    virtual void ReadElement(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        Element::Pointer& pThisElement
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadElement). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads an array of elements
     * @param rThisNodes The nodes constituying the element
     * @param rThisProperties The Properties of the element
     * @param rThisElement The array of elements
     */
    virtual void ReadElements(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        ElementsContainerType& rThisElements
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadElements). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the elements connectivities
     * @param rElementsConnectivities The elements connectivities
     * @return The number of elements
     */
    virtual std::size_t ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
    {
        KRATOS_ERROR << "Calling base class method (ReadElementsConnectivities). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes an array of elements
     * @param rThisElements The array of elements to be written
     */
    virtual void WriteElements(ElementsContainerType const& rThisElements)
    {
        KRATOS_ERROR << "Calling base class method (WriteElements). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads one condition
     * @param rThisNodes The nodes constituying the condition
     * @param rThisProperties The Properties of the condition
     * @param pThisCondition The pointer to the condition
     */
    virtual void ReadCondition(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        Condition::Pointer& pThisCondition
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadCondition). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads an array of conditions
     * @param rThisNodes The nodes constituying the condition
     * @param rThisProperties The Properties of the condition
     * @param rThisConditions The array of conditions
     */
    virtual void ReadConditions(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        ConditionsContainerType& rThisConditions
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadConditions). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the conditions connectivities
     * @param rConditionsConnectivities The conditions connectivities
     * @return The number of conditions
     */
    virtual std::size_t ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
    {
        KRATOS_ERROR << "Calling base class method (ReadConditionsConnectivities). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes an array of conditions
     * @param rThisConditions The array of conditions to be written
     */
    virtual void WriteConditions(ConditionsContainerType const& rThisConditions)
    {
        KRATOS_ERROR << "Calling base class method (WriteConditions). Please check the definition of derived class" << std::endl;
    }
    
    /**
     * @brief Reads the master-slave constraints from an input source.
     * @details This method is intended to be overridden by derived classes to implement
     * the specific logic for reading master-slave constraints into the provided
     * container. The base class implementation throws an error, indicating that
     * the method must be implemented in the derived class.
     * @param rThisNodes The nodes to be used for associating the master-slave constraints.
     * @param rConstraintContainer The container where the master-slave
     *        constraints will be stored. This container is expected to be populated
     *        by the derived class implementation.
     * @throws Exception If the base class method is called directly, an error is
     *         thrown to indicate that the method must be implemented in a derived class.
     */
    virtual void ReadConstraints(
        NodesContainerType& rThisNodes,
        MasterSlaveConstraintContainerType& rConstraintContainer
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadConstraints). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the constraints connectivities
     * @param rConditionsConnectivities The constraints connectivities
     * @return The number of constraints
     */
    virtual std::size_t ReadConstraintsConnectivities(ConnectivitiesContainerType& rConstraintsConnectivities)
    {
        KRATOS_ERROR << "Calling base class method (ReadConstraintsConnectivities). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief Writes the master-slave constraints to the output.
     * @details This method is intended to be overridden by derived classes to provide
     * specific functionality for writing master-slave constraints. The base
     * class implementation throws an error, indicating that the method must
     * be implemented in the derived class.
     * @param rConstraintContainer The container holding the master-slave
     *        constraints to be written.
     * @throws Exception Always throws an error if called on the base class.
     *         Derived classes must override this method to provide the actual
     *         implementation.
     */
    virtual void WriteConstraints(MasterSlaveConstraintContainerType const& rConstraintContainer)
    {
        KRATOS_ERROR << "Calling base class method (WriteNewConstraint). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the initial values of the model part
     * @param rThisModelPart The model part with the initial values to be read
     */
    virtual void ReadInitialValues(ModelPart& rThisModelPart)
    {
        KRATOS_ERROR << "Calling base class method (ReadInitialValues). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the initial values of the nodes, elements and conditios
     * @param rThisNodes The nodes with the initial values to be read
     * @param rThisElements The elements with the initial values to be read
     * @param rThisConditions The conditions with the initial values to be read
     */
    virtual void ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions)
    {
        KRATOS_ERROR << "Calling base class method (ReadInitialValues). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method reads the mesh
     * @param rThisMesh The mesh to be read
     */
    virtual void ReadMesh(MeshType & rThisMesh)
    {
        KRATOS_ERROR << "Calling base class method (ReadMesh). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the mesh
     * @param rThisMesh The mesh to be written
     */
    KRATOS_DEPRECATED_MESSAGE("'WriteMesh' with a non-const Mesh as input is deprecated. Please use the version of this function that accepts a const Mesh instead.")
    virtual void WriteMesh( MeshType& rThisMesh )
    {
        KRATOS_ERROR << "Calling base class method (WriteMesh). Please check the implementation of derived classes" << std::endl;
    }

    /**
     * @brief This method writes the mesh
     * @param rThisMesh The mesh to be written
     */
    virtual void WriteMesh(const MeshType& rThisMesh )
    {
        // legacy for backward compatibility
        MeshType& non_const_mesh = const_cast<MeshType&>(rThisMesh);
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->WriteMesh(non_const_mesh);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING

        // activate this error once the legacy code is removed
        // KRATOS_ERROR << "Calling base class method (WriteMesh). Please check the implementation of derived classes" << std::endl;
    }

    /**
     * @brief This method reads the model part
     * @param rThisModelPart The model part to be read
     */
    virtual void ReadModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_ERROR << "Calling base class method (ReadModelPart). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the model part
     * @param rThisModelPart The model part to be written
     */
    KRATOS_DEPRECATED_MESSAGE("'WriteModelPart' with a non-const ModelPart as input is deprecated. Please use the version of this function that accepts a const ModelPart instead.")
    virtual void WriteModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_ERROR << "Calling base class method (WriteModelPart). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the model part
     * @param rThisModelPart The model part to be written
     */
    virtual void WriteModelPart(const ModelPart& rThisModelPart)
    {
        // legacy for backward compatibility
        ModelPart& non_const_model_part = const_cast<ModelPart&>(rThisModelPart);
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->WriteModelPart(non_const_model_part);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING

        // activate this error once the legacy code is removed
        // KRATOS_ERROR << "Calling base class method (WriteModelPart). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the node mesh
     * @param rThisMesh The mesh to be written
     */
    KRATOS_DEPRECATED_MESSAGE("'WriteNodeMesh' with a non-const Mesh as input is deprecated. Please use the version of this function that accepts a const Mesh instead.")
    virtual void WriteNodeMesh( MeshType& rThisMesh )
    {
        KRATOS_ERROR << "Calling base class method (WriteNodeMesh). Please check the implementation of derived classes" << std::endl;
    }

    /**
     * @brief This method writes the node mesh
     * @param rThisMesh The mesh to be written
     */
    virtual void WriteNodeMesh(const MeshType& rThisMesh )
    {
        // legacy for backward compatibility
        MeshType& non_const_mesh = const_cast<MeshType&>(rThisMesh);
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->WriteNodeMesh(non_const_mesh);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING

        // activate this error once the legacy code is removed
        // KRATOS_ERROR << "Calling base class method (WriteNodeMesh). Please check the implementation of derived classes" << std::endl;
    }

    /**
     * @brief Read the input file and create the nodal connectivities graph, stored in CSR format.
     * @details This function produces input for Metis' nodal graph partitioning algorithms.
     * The nodal graph is stored as a (compressed) matrix where index (i,j) is non-zero if
     * there is an edge in the mesh joining nodes i and j (note that nodes are numbered from zero here,
     * to make integration with Metis simpler).
     * After call, will point to C array of size NumNodes+1 containing the first CSR array: entries related to node k are stored between positions (*NodeIndices)[k] and (*NodeIndices)[k+1] of *NodeConnectivities.
     * @param rAuxConnectivities After call, will point to a C array of size (*NodeIndices)[NumNodes].
     * entries between (*NodeIndices)[k] and (*NodeIndices)[k+1] are a list of all nodes connected
     * to node k (counting from 0).
     * @return Number of nodes.
     */
    virtual std::size_t ReadNodalGraph(ConnectivitiesContainerType& rAuxConnectivities)
    {
        KRATOS_ERROR << "Calling base class method (ReadNodalGraph). Please check the definition of derived class" << std::endl;;
    }

    /**
     * @brief This method divides a model part into partitions
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    virtual void DivideInputToPartitions(SizeType NumberOfPartitions,
                                         const PartitioningInfo& rPartitioningInfo)
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        DivideInputToPartitions(NumberOfPartitions, rPartitioningInfo.Graph, rPartitioningInfo.NodesPartitions, rPartitioningInfo.ElementsPartitions, rPartitioningInfo.ConditionsPartitions, rPartitioningInfo.NodesAllPartitions, rPartitioningInfo.ElementsAllPartitions, rPartitioningInfo.ConditionsAllPartitions); // for backward compatibility
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING

        // KRATOS_ERROR << "Calling base class method (DivideInputToPartitions). Please check the definition of derived class" << std::endl; // enable this once the old version of this function is removed
    }

    /**
     * @brief This method divides a model part into partitions
     * @param NumberOfPartitions The number of partitions
     * @param rDomainsColoredGraph The colors of the partition graph
     * @param rNodesPartitions The partitions indices of the nodes
     * @param rElementsPartitions The partitions indices of the elements
     * @param rConditionsPartitions The partitions indices of the conditions
     * @param rNodesAllPartitions The partitions of the nodes
     * @param rElementsAllPartitions The partitions of the elements
     * @param rConditionsAllPartitions The partitions of the conditions
     */
    KRATOS_DEPRECATED_MESSAGE("'This version of \"DivideInputToPartitions\" is deprecated, please use the interface that accepts a \"PartitioningInfo\"")
    virtual void DivideInputToPartitions(SizeType NumberOfPartitions,
                                         GraphType const& rDomainsColoredGraph,
                                         PartitionIndicesType const& rNodesPartitions,
                                         PartitionIndicesType const& rElementsPartitions,
                                         PartitionIndicesType const& rConditionsPartitions,
                                         PartitionIndicesContainerType const& rNodesAllPartitions,
                                         PartitionIndicesContainerType const& rElementsAllPartitions,
                                         PartitionIndicesContainerType const& rConditionsAllPartitions)
    {
        KRATOS_ERROR << "Calling base class method (DivideInputToPartitions). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method divides a model part into partitions
     * @param pStreams The stream pointer
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    virtual void DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                         SizeType NumberOfPartitions,
                                         const PartitioningInfo& rPartitioningInfo)
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        DivideInputToPartitions(pStreams, NumberOfPartitions, rPartitioningInfo.Graph, rPartitioningInfo.NodesPartitions, rPartitioningInfo.ElementsPartitions, rPartitioningInfo.ConditionsPartitions, rPartitioningInfo.NodesAllPartitions, rPartitioningInfo.ElementsAllPartitions, rPartitioningInfo.ConditionsAllPartitions); // for backward compatibility
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING

        // KRATOS_ERROR << "Calling base class method (DivideInputToPartitions). Please check the definition of derived class" << std::endl; // enable this once the old version of this function is removed
    }

    /**
     * @brief This method divides a model part into partitions
     * @param pStreams The stream pointer
     * @param NumberOfPartitions The number of partitions
     * @param rDomainsColoredGraph The colors of the partition graph
     * @param rNodesPartitions The partitions indices of the nodes
     * @param rElementsPartitions The partitions indices of the elements
     * @param rConditionsPartitions The partitions indices of the conditions
     * @param rNodesAllPartitions The partitions of the nodes
     * @param rElementsAllPartitions The partitions of the elements
     * @param rConditionsAllPartitions The partitions of the conditions
     */
    KRATOS_DEPRECATED_MESSAGE("'This version of \"DivideInputToPartitions\" is deprecated, please use the interface that accepts a \"PartitioningInfo\"")
    virtual void DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                         SizeType NumberOfPartitions,
                                         GraphType const& rDomainsColoredGraph,
                                         PartitionIndicesType const& rNodesPartitions,
                                         PartitionIndicesType const& rElementsPartitions,
                                         PartitionIndicesType const& rConditionsPartitions,
                                         PartitionIndicesContainerType const& rNodesAllPartitions,
                                         PartitionIndicesContainerType const& rElementsAllPartitions,
                                         PartitionIndicesContainerType const& rConditionsAllPartitions)
    {
        KRATOS_ERROR << "Calling base class method (DivideInputToPartitions). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief Virtual method to read element and condition IDs from a sub-model part.
     * @details This method is intended to be overridden by derived classes.
     * @param rModelPartName The name of the sub-model part to read from.
     * @param rElementsIds Set to store element IDs.
     * @param rConditionsIds Set to store condition IDs.
     */
    virtual void ReadSubModelPartElementsAndConditionsIds(
        std::string const& rModelPartName,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadSubModelPartElementsAndConditionsIds). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief Virtual method to read element, condition, master-slave constraint, and geometry IDs from a sub-model part.
     * @details This method is intended to be overridden by derived classes.
     * @param rModelPartName The name of the sub-model part to read from.
     * @param rElementsIds Set to store element IDs.
     * @param rConditionsIds Set to store condition IDs.
     * @param rConstraintIds Set to store master-slave constraint IDs.
     * @param rGeometriesIds Set to store geometry IDs.
     */
    virtual void ReadSubModelPartEntitiesIds(
        std::string const& rModelPartName,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds,
        std::unordered_set<SizeType>& rConstraintIds,
        std::unordered_set<SizeType>& rGeometriesIds
        )
    {
        KRATOS_WARNING("IO") << " The method ReadSubModelPartEntitiesIds with Constraint and Geometries is not implemented. Only the elements and conditions are read." << std::endl;
        return ReadSubModelPartElementsAndConditionsIds(rModelPartName, rElementsIds, rConditionsIds);
    }

    /**
     * @brief Virtual method to read nodal graph from entities list.
     * @details This method is intended to be overridden by derived classes.
     * @param rAuxConnectivities Container of connectivities.
     * @param rElementsIds Set of element IDs.
     * @param rConditionsIds Set of condition IDs.
     * @return The size of the nodal graph.
     */
    virtual std::size_t ReadNodalGraphFromEntitiesList(
        ConnectivitiesContainerType& rAuxConnectivities,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds
        )
    {
        KRATOS_ERROR << "Calling base class method (ReadNodalGraphFromEntitiesList). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief Virtual method to read nodal graph from entities list including master-slave constraints and geometries.
     * @details This method is intended to be overridden by derived classes.
     * @param rAuxConnectivities Container of connectivities.
     * @param rElementsIds Set of element IDs.
     * @param rConditionsIds Set of condition IDs.
     * @param rConstraintIds Set of master-slave constraint IDs.
     * @param rGeometriesIds Set of geometry IDs.
     * @return The size of the nodal graph.
     */
    virtual std::size_t ReadNodalGraphFromEntitiesList(
        ConnectivitiesContainerType& rAuxConnectivities,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds,
        std::unordered_set<SizeType>& rConstraintIds,
        std::unordered_set<SizeType>& rGeometriesIds
        )
    {
        KRATOS_WARNING("IO") << " The method ReadNodalGraphFromEntitiesList with Constraint and Geometries is not implemented. Only the elements and conditions are read." << std::endl;
        return ReadNodalGraphFromEntitiesList(rAuxConnectivities, rElementsIds, rConditionsIds);
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
        return "IO";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IO";
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

}; // Class IO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                IO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const IO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.
