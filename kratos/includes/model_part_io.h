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
//                   Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <filesystem>
#include <fstream>

// External includes

// Project includes
#include "includes/io.h"

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

/// An IO class for reading and writing a modelpart
/** This class reads and writes all modelpart data including the meshes.
*/
class KRATOS_API(KRATOS_CORE) ModelPartIO 
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPartIO
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartIO);

    /// Alias for the base IO type.
    using BaseType = IO;

    /// Alias for the node type.
    using NodeType = BaseType::NodeType;

    /// Alias for the mesh type.
    using MeshType = BaseType::MeshType;

    /// Alias for the nodes container type.
    using NodesContainerType = BaseType::NodesContainerType;

    /// Alias for the properties container type.
    using PropertiesContainerType = BaseType::PropertiesContainerType;

    /// Alias for the elements container type.
    using ElementsContainerType = BaseType::ElementsContainerType;

    /// Alias for the conditions container type.
    using ConditionsContainerType = BaseType::ConditionsContainerType;

    /// Alias for the master-slave constraints container type.
    using MasterSlaveConstraintContainerType = BaseType::MasterSlaveConstraintContainerType;

    /// Alias for the connectivities container type.
    using ConnectivitiesContainerType = BaseType::ConnectivitiesContainerType;

    /// Alias for the output files container type.
    using OutputFilesContainerType = std::vector<std::ostream*>;

    /// Alias for the size type.
    using SizeType = std::size_t;

    /// Prevents this class from hiding IO::WriteProperties(Properties)
    using BaseType::WriteProperties;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with filename.
    ModelPartIO(
        std::filesystem::path const& Filename,
        const Flags Options = IO::READ | IO::IGNORE_VARIABLES_ERROR.AsFalse() | IO::SKIP_TIMER);

    /// Constructor with stream.
    ModelPartIO(
        Kratos::shared_ptr<std::iostream> Stream,
        const Flags Options = IO::IGNORE_VARIABLES_ERROR.AsFalse() | IO::SKIP_TIMER);

    /// Constructor with filenames.
    // ModelPartIO(std::string const& InputFilename, std::string const& OutputFilename)
    //     : mNumberOfLines(0), mInput(std::ifstream(InputFilename.c_str())), mOutput(std::ofstream(OutputFilename.c_str()))
    // {
    // }

    /// Destructor.
    ~ModelPartIO() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method reads one node
     * @param rThisNode The node to be read
     */
    bool ReadNode(NodeType& rThisNode) override;

    /**
     * @brief This method reads the nodes from an array of nodes
     * @param rThisNodes The array of nodes to be read
     */
    bool ReadNodes(NodesContainerType& rThisNodes) override;

    /**
     * @brief This method reads the number of nodes
     * @return The number of nodes
     */
    std::size_t ReadNodesNumber() override;

    /**
     * @brief This method writes the nodes from an array of nodes
     * @param rThisNodes The array of nodes to be written
     */
    void WriteNodes(NodesContainerType const& rThisNodes) override;

    /**
     * @brief This method reads one Properties
     * @param rThisProperties The Properties to be read
     */
    void ReadProperties(Properties& rThisProperties) override;

    /**
     * @brief This method reads the Properties from an array of Properties
     * @param rThisProperties The array of Properties to be read
     */
    void ReadProperties(PropertiesContainerType& rThisProperties) override;

    /**
     * @brief This method writes one Properties
     * @param rThisProperties The Properties to be written
     */
    void WriteProperties(PropertiesContainerType const& rThisProperties) override;

    /**
     * @brief This method reads one geometry
     * @param rThisNodes The nodes constituying the geometry
     * @param pThisGeometries The pointer to the geometry
     */
    void ReadGeometry(
        NodesContainerType& rThisNodes,
        GeometryType::Pointer& pThisGeometry
        ) override;

    /**
     * @brief This method reads an array of geometries
     * @param rThisNodes The nodes constituying the geometry
     * @param rThisGeometry The array of geometries
     */
    void ReadGeometries(
        NodesContainerType& rThisNodes,
        GeometryContainerType& rThisGeometries
        ) override;

    /**
     * @brief This method reads the geometries connectivities
     * @param rGeometriesConnectivities The geometries connectivities
     * @return The number of geometries
     */
    std::size_t ReadGeometriesConnectivities(ConnectivitiesContainerType& rGeometriesConnectivities) override;

    /**
     * @brief This method writes an array of geometries
     * @param rThisGeometries The array of geometries to be written
     */
    void WriteGeometries(GeometryContainerType const& rThisGeometries) override;

    /**
     * @brief This method reads one element
     * @param rThisNodes The nodes constituying the element
     * @param rThisProperties The Properties of the element
     * @param pThisElements The pointer to the element
     */
    void ReadElement(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        Element::Pointer& pThisElement
        ) override;

    /**
     * @brief This method reads an array of elements
     * @param rThisNodes The nodes constituying the element
     * @param rThisProperties The Properties of the element
     * @param rThisElement The array of elements
     */
    void ReadElements(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        ElementsContainerType& rThisElements
        ) override;

    /**
     * @brief This method reads the elements connectivities
     * @param rElementsConnectivities The elements connectivities
     * @return The number of elements
     */
    std::size_t ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities) override;

    /**
     * @brief This method writes an array of elements
     * @param rThisElements The array of elements to be written
     */
    void WriteElements(ElementsContainerType const& rThisElements) override;

    /**
     * @brief This method reads an array of conditions
     * @param rThisNodes The nodes constituying the condition
     * @param rThisProperties The Properties of the condition
     * @param rThisConditions The array of conditions
     */
    void ReadConditions(
        NodesContainerType& rThisNodes,
        PropertiesContainerType& rThisProperties,
        ConditionsContainerType& rThisConditions
        ) override;

    /**
     * @brief This method reads the conditions connectivities
     * @param rConditionsConnectivities The conditions connectivities
     * @return The number of conditions
     */
    std::size_t ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities) override;

    /**
     * @brief This method writes an array of conditions
     * @param rThisConditions The array of conditions to be written
     */
    void WriteConditions(ConditionsContainerType const& rThisConditions) override;

    /**
     * @brief Reads the master-slave constraints from an input source.
     * @details This method is intended to be overridden by derived classes to implement
     * the specific logic for reading master-slave constraints into the provided
     * container. The base class implementation throws an error, indicating that
     * the method must be implemented in the derived class.
     * @param rThisNodes The nodes to be used for associating the master-slave constraints.
     * @param rMasterSlaveConstraintContainer The container where the master-slave
     *        constraints will be stored. This container is expected to be populated
     *        by the derived class implementation.
     */
    void ReadMasterSlaveConstraints(
        NodesContainerType& rThisNodes,
        MasterSlaveConstraintContainerType& rMasterSlaveConstraintContainer
        ) override;

    /**
     * @brief Writes the master-slave constraints to the output.
     * @details This method is intended to be overridden by derived classes to provide
     * specific functionality for writing master-slave constraints. The base
     * class implementation throws an error, indicating that the method must
     * be implemented in the derived class.
     * @param rMasterSlaveConstraintContainer The container holding the master-slave
     *        constraints to be written.
     */
    void WriteMasterSlaveConstraints(MasterSlaveConstraintContainerType const& rMasterSlaveConstraintContainer) override;

    /**
     * @brief This method reads the initial values of the model part
     * @param rThisModelPart The model part with the initial values to be read
     */
    void ReadInitialValues(ModelPart& rThisModelPart) override;

    /**
     * @brief This method reads the mesh
     * @param rThisMesh The mesh to be read
     */
    void ReadMesh(MeshType & rThisMesh) override;

    /**
     * @brief This method writes the mesh
     * @param rThisMesh The mesh to be written
     */
    void WriteMesh(MeshType & rThisMesh) override;

    /**
     * @brief This method reads the model part
     * @param rThisModelPart The model part to be read
     */
    void ReadModelPart(ModelPart & rThisModelPart) override;

    /**
     * @brief This method writes the model part
     * @param rThisModelPart The model part to be written
     */
    void WriteModelPart(ModelPart & rThisModelPart) override;

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
    std::size_t ReadNodalGraph(ConnectivitiesContainerType& rAuxConnectivities) override;

    /**
     * @brief This method divides a model part into partitions
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    void DivideInputToPartitions(SizeType NumberOfPartitions,
                                 const PartitioningInfo& rPartitioningInfo) override;

    /**
     * @brief This method divides a model part into partitions
     * @param pStreams The stream pointer
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    void DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                SizeType NumberOfPartitions,
                                const PartitioningInfo& rPartitioningInfo) override;

    void SwapStreamSource(Kratos::shared_ptr<std::iostream> newStream);

    /**
     * @brief Virtual method to read element and condition IDs from a sub-model part.
     * @details This method is intended to be overridden by derived classes.
     * @param rModelPartName The name of the sub-model part to read from.
     * @param rElementsIds Set to store element IDs.
     * @param rConditionsIds Set to store condition IDs.
     */
    void ReadSubModelPartElementsAndConditionsIds(
        std::string const& rModelPartName,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds
        ) override;

    /**
     * @brief Virtual method to read element, condition, master-slave constraint, and geometry IDs from a sub-model part.
     * @details This method is intended to be overridden by derived classes.
     * @param rModelPartName The name of the sub-model part to read from.
     * @param rElementsIds Set to store element IDs.
     * @param rConditionsIds Set to store condition IDs.
     * @param rMasterSlaveConstraintIds Set to store master-slave constraint IDs.
     * @param rGeometriesIds Set to store geometry IDs.
     */
    void ReadSubModelPartElementsAndConditionsIds(
        std::string const& rModelPartName,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds,
        std::unordered_set<SizeType>& rMasterSlaveConstraintIds,
        std::unordered_set<SizeType>& rGeometriesIds
        ) override;

    /**
     * @brief Virtual method to read nodal graph from entities list.
     * @details This method is intended to be overridden by derived classes.
     * @param rAuxConnectivities Container of connectivities.
     * @param rElementsIds Set of element IDs.
     * @param rConditionsIds Set of condition IDs.
     * @return The size of the nodal graph.
     */
    std::size_t ReadNodalGraphFromEntitiesList(
        ConnectivitiesContainerType& rAuxConnectivities,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds
        ) override;

    /**
     * @brief Virtual method to read nodal graph from entities list including master-slave constraints and geometries.
     * @details This method is intended to be overridden by derived classes.
     * @param rAuxConnectivities Container of connectivities.
     * @param rElementsIds Set of element IDs.
     * @param rConditionsIds Set of condition IDs.
     * @param rMasterSlaveConstraintIds Set of master-slave constraint IDs.
     * @param rGeometriesIds Set of geometry IDs.
     * @return The size of the nodal graph.
     */
    std::size_t ReadNodalGraphFromEntitiesList(
        ConnectivitiesContainerType& rAuxConnectivities,
        std::unordered_set<SizeType>& rElementsIds,
        std::unordered_set<SizeType>& rConditionsIds,
        std::unordered_set<SizeType>& rMasterSlaveConstraintIds,
        std::unordered_set<SizeType>& rGeometriesIds
        ) override;

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
        return "ModelPartIO";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ModelPartIO";
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

    /**
     * @brief Returns the reordered node ID.
     * @param NodeId The original node ID.
     * @return The reordered node ID (same as input in this base class).
     */
    virtual ModelPartIO::SizeType ReorderedNodeId(ModelPartIO::SizeType NodeId);

    /**
     * @brief Returns the reordered geometry ID.
     * @param GeometryId The original geometry ID.
     * @return The reordered geometry ID (same as input in this base class).
     */
    virtual ModelPartIO::SizeType ReorderedGeometryId(ModelPartIO::SizeType GeometryId);

    /**
     * @brief Returns the reordered element ID.
     * @param ElementId The original element ID.
     * @return The reordered element ID (same as input in this base class).
     */
    virtual ModelPartIO::SizeType ReorderedElementId(ModelPartIO::SizeType ElementId);

    /**
     * @brief Returns the reordered condition ID.
     * @param ConditionId The original condition ID.
     * @return The reordered condition ID (same as input in this base class).
     */
    virtual ModelPartIO::SizeType ReorderedConditionId(ModelPartIO::SizeType ConditionId);

    /**
     * @brief Returns the reordered master-slave constraint ID.
     * @param ConstraintId The original constraint ID.
     * @return The reordered constraint ID (same as input in this base class).
     */
    virtual ModelPartIO::SizeType ReorderedMasterSlaveConstraintId(ModelPartIO::SizeType ConstraintId);

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

    SizeType mNumberOfLines;

    std::filesystem::path mBaseFilename;
    Flags mOptions;

    Kratos::shared_ptr<std::iostream> mpStream;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    std::string& ReadBlockName(std::string& rBlockName);

    void SkipBlock(std::string const& BlockName);

    bool CheckEndBlock(std::string const& BlockName, std::string& rWord);

    void ReadModelPartDataBlock(ModelPart& rModelPart, const bool IsSubmodelpart=false);

    void WriteModelPartDataBlock(ModelPart& rModelPart, const bool IsSubmodelpart=false);

    template<class TablesContainerType>
    void ReadTableBlock(TablesContainerType& rTables);

    void ReadTableBlock(ModelPart::TablesContainerType& rTables);

    template<class TablesContainerType>
    void WriteTableBlock(TablesContainerType& rTables);

    void WriteTableBlock(ModelPart::TablesContainerType& rTables);

    void ReadNodesBlock(NodesContainerType& rThisNodes);

    void ReadNodesBlock(ModelPart& rModelPart);

    std::size_t CountNodesInBlock();

    void ReadPropertiesBlock(PropertiesContainerType& rThisProperties);

    void ReadGeometriesBlock(ModelPart& rModelPart);

    void ReadGeometriesBlock(NodesContainerType& rThisNodes, GeometryContainerType& rThisGeometries);

    void ReadElementsBlock(ModelPart& rModelPart);

    void ReadElementsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements);

    void ReadConditionsBlock(ModelPart& rModelPart);

    void ReadConditionsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions);

    /**
     * @brief Reads the master-slave constraints block from the input file.
     * @details This function is responsible for parsing and loading the master-slave 
     * constraints defined in the input file into the provided ModelPart. 
     * Master-slave constraints are used to define relationships between 
     * degrees of freedom in the model.
     * @param rModelPart Reference to the ModelPart where the constraints will be added.
     */
    void ReadMasterSlaveConstraintsBlock(ModelPart& rModelPart);

    /**
     * @brief Reads a block of master-slave constraints from the input.
     * @details This function processes and reads the master-slave constraint data
     * from the input file, associating it with the provided nodes and
     * storing the constraints in the specified container.
     * @param rThisNodes A reference to the container of nodes to be used
     *                   for associating the master-slave constraints.
     * @param rMasterSlaveConstraints A reference to the container where
     *                                the read master-slave constraints
     *                                will be stored.
     */
    void ReadMasterSlaveConstraintsBlock(
        NodesContainerType& rThisNodes,
        MasterSlaveConstraintContainerType& rMasterSlaveConstraints
        );

    /**
     * @brief Reads nodal data block from the input stream into the given model part.
     * @param rThisModelPart Reference to the model part whose nodal data will be read.
     */
    void ReadNodalDataBlock(ModelPart& rThisModelPart);

    /**
     * @brief Writes nodal data block from the given model part to the output stream.
     * @param rThisModelPart Reference to the model part whose nodal data will be written.
     */
    void WriteNodalDataBlock(ModelPart& rThisModelPart);

    /**
     * @brief Reads nodal DOF variable data from the input stream into the provided nodes.
     * @tparam TVariableType Type of the variable to be read.
     * @param rThisNodes Reference to the container of nodes.
     * @param rVariable The variable to be read.
     */
    template<class TVariableType>
    void ReadNodalDofVariableData(NodesContainerType& rThisNodes, const TVariableType& rVariable);

    /**
     * @brief Reads nodal flags from the input stream and assigns them to the provided nodes.
     * @param rThisNodes Reference to the container of nodes.
     * @param rFlags Flags to be read and applied.
     */
    void ReadNodalFlags(NodesContainerType& rThisNodes, Flags const& rFlags);

    /**
     * @brief Reads scalar variable data for each node in the container.
     * @tparam TVariableType Type of the scalar variable.
     * @param rThisNodes Reference to the container of nodes.
     * @param rVariable The scalar variable to read.
     */
    template<class TVariableType>
    void ReadNodalScalarVariableData(NodesContainerType& rThisNodes, const TVariableType& rVariable);

    /**
     * @brief Reads vectorial variable data for each node in the container.
     * @tparam TVariableType Type of the variable.
     * @tparam TDataType Type of the vector data.
     * @param rThisNodes Reference to the container of nodes.
     * @param rVariable The variable to read.
     * @param Dummy Dummy parameter to help with template deduction.
     */
    template<class TVariableType, class TDataType>
    void ReadNodalVectorialVariableData(NodesContainerType& rThisNodes, const TVariableType& rVariable, TDataType Dummy);

    /**
     * @brief Reads elemental data block from the input stream into the given elements container.
     * @param rThisElements Reference to the container of elements.
     */
    void ReadElementalDataBlock(ElementsContainerType& rThisElements);

    /**
     * @brief Reads scalar variable data for each element in the container.
     * @tparam TVariableType Type of the scalar variable.
     * @param rThisElements Reference to the container of elements.
     * @param rVariable The scalar variable to read.
     */
    template<class TVariableType>
    void ReadElementalScalarVariableData(ElementsContainerType& rThisElements, const TVariableType& rVariable);

    /**
     * @brief Reads vectorial variable data for each element in the container.
     * @tparam TVariableType Type of the variable.
     * @tparam TDataType Type of the vector data.
     * @param rThisElements Reference to the container of elements.
     * @param rVariable The variable to read.
     * @param Dummy Dummy parameter to help with template deduction.
     */
    template<class TVariableType, class TDataType>
    void ReadElementalVectorialVariableData(ElementsContainerType& rThisElements, const TVariableType& rVariable, TDataType Dummy);

    /**
     * @brief Reads conditional data block from the input stream into the given conditions container.
     * @param rThisConditions Reference to the container of conditions.
     */
    void ReadConditionalDataBlock(ConditionsContainerType& rThisConditions);

    /**
     * @brief Reads scalar variable data for each condition in the container.
     * @tparam TVariableType Type of the scalar variable.
     * @param rThisConditions Reference to the container of conditions.
     * @param rVariable The scalar variable to read.
     */
    template<class TVariableType>
    void ReadConditionalScalarVariableData(ConditionsContainerType& rThisConditions, const TVariableType& rVariable);

    /**
     * @brief Reads vectorial variable data for each condition in the container.
     * @tparam TVariableType Type of the variable.
     * @tparam TDataType Type of the vector data.
     * @param rThisConditions Reference to the container of conditions.
     * @param rVariable The variable to read.
     * @param Dummy Dummy parameter to help with template deduction.
     */
    template<class TVariableType, class TDataType>
    void ReadConditionalVectorialVariableData(ConditionsContainerType& rThisConditions, const TVariableType& rVariable, TDataType Dummy);

    /**
     * @brief Reads master-slave constraint data block from the input stream into the given constraints container.
     * @param rThisConstraints Reference to the container of master-slave constraints.
     */
    void ReadMasterSlaveConstraintDataBlock(MasterSlaveConstraintContainerType& rThisConstraints);

    /**
     * @brief Reads scalar variable data for each constraint in the container.
     * @tparam TVariableType Type of the scalar variable.
     * @param rThisConstraints Reference to the container of constraints.
     * @param rVariable The scalar variable to read.
     */
    template<class TVariableType>
    void ReadMasterSlaveConstraintScalarVariableData(MasterSlaveConstraintContainerType& rThisConstraints, const TVariableType& rVariable);

    /**
     * @brief Reads vectorial variable data for each constraint in the container.
     * @tparam TVariableType Type of the variable.
     * @tparam TDataType Type of the vector data.
     * @param rThisConstraints Reference to the container of constraints.
     * @param rVariable The variable to read.
     * @param Dummy Dummy parameter to help with template deduction.
     */
    template<class TVariableType, class TDataType>
    void ReadMasterSlaveConstraintVectorialVariableData(MasterSlaveConstraintContainerType& rThisConstraints, const TVariableType& rVariable, TDataType Dummy);

    /**
     * @brief Writes a data block for the specified object container.
     * @tparam TObjectsContainerType Type of the object container.
     * @param rThisObjectContainer Reference to the container of objects.
     * @param rObjectName Name of the object group to write.
     */
    template<class TObjectsContainerType>
    void WriteDataBlock(const TObjectsContainerType& rThisObjectContainer, const std::string& rObjectName);

    /**
     * @brief Writes a data block for the specified variable in the object container.
     * @tparam TVariableType Type of the variable.
     * @tparam TObjectsContainerType Type of the object container.
     * @param rThisObjectContainer Reference to the container of objects.
     * @param rVariable Pointer to the variable data to write.
     * @param rObjectName Name of the object group to write.
     */
    template<class TVariableType, class TObjectsContainerType>
    void WriteDataBlock(const TObjectsContainerType& rThisObjectContainer, const VariableData* rVariable, const std::string& rObjectName);

    /**
     * @brief Reads the geometries connectivities block from the input stream.
     * @param rThisConnectivities Reference to the container where the connectivities will be stored.
     * @return The number of geometries read.
     */
    SizeType ReadGeometriesConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);

    /**
     * @brief Reads the elements connectivities block from the input stream.
     * @param rThisConnectivities Reference to the container where the connectivities will be stored.
     * @return The number of elements read.
     */
    SizeType ReadElementsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);

    /**
     * @brief Reads the conditions connectivities block from the input stream.
     * @param rThisConnectivities Reference to the container where the connectivities will be stored.
     * @return The number of conditions read.
     */
    SizeType ReadConditionsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);

    /**
     * @brief Reads the master-slave constraints connectivities block from the input stream.
     * @param rThisConnectivities Reference to the container where the connectivities will be stored.
     * @return The number of master-slave constraints read.
     */
    SizeType ReadMasterSlaveConstraintsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);

    /**
     * @brief Fills the nodal connectivities from the geometry block.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     */
    void FillNodalConnectivitiesFromGeometryBlock(ConnectivitiesContainerType& rNodalConnectivities);

    /**
     * @brief Fills the nodal connectivities from the element block.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     */
    void FillNodalConnectivitiesFromElementBlock(ConnectivitiesContainerType& rNodalConnectivities);

    /**
     * @brief Fills the nodal connectivities from the condition block.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     */
    void FillNodalConnectivitiesFromConditionBlock(ConnectivitiesContainerType& rNodalConnectivities);

    /**
     * @brief Fills the nodal connectivities from the master-slave constraint block.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     */
    void FillNodalConnectivitiesFromMasterSlaveConstraintBlock(ConnectivitiesContainerType& rNodalConnectivities);

    /**
     * @brief Fills the nodal connectivities from the geometry block for a specific list of geometries.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     * @param rGeometriesIds Set of geometry IDs to filter the connectivities.
     */
    void FillNodalConnectivitiesFromGeometryBlockInList(
        ConnectivitiesContainerType& rNodalConnectivities,
        std::unordered_set<SizeType>& rGeometriesIds
        );

    /**
     * @brief Fills the nodal connectivities from the element block for a specific list of elements.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     * @param rElementsIds Set of element IDs to filter the connectivities.
     */
    void FillNodalConnectivitiesFromElementBlockInList(
        ConnectivitiesContainerType& rNodalConnectivities,
        std::unordered_set<SizeType>& rElementsIds
        );

    /**
     * @brief Fills the nodal connectivities from the condition block for a specific list of conditions.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     * @param rConditionsIds Set of condition IDs to filter the connectivities.
     */
    void FillNodalConnectivitiesFromConditionBlockInList(
        ConnectivitiesContainerType& rNodalConnectivities,
        std::unordered_set<SizeType>& rConditionsIds
        );

    /**
     * @brief Fills the nodal connectivities from the master-slave constraint block for a specific list of constraints.
     * @param rNodalConnectivities Reference to the container where the nodal connectivities will be stored.
     * @param rConstraintsIds Set of master-slave constraint IDs to filter the connectivities.
     */
    void FillNodalConnectivitiesFromMasterSlaveConstraintBlockInList(
        ConnectivitiesContainerType& rNodalConnectivities,
        std::unordered_set<SizeType>& rConstraintsIds
        );

    void ReadCommunicatorDataBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);

    void ReadCommunicatorLocalNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);

    void ReadCommunicatorGhostNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);

    void ReadMeshBlock(ModelPart& rModelPart);

    void WriteMeshBlock(ModelPart& rModelPart);

    void ReadMeshDataBlock(MeshType& rMesh);

    void ReadMeshNodesBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadMeshElementsBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadMeshConditionsBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadMeshPropertiesBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadSubModelPartBlock(ModelPart& rMainModelPart, ModelPart& rParentModelPart);

    void WriteSubModelPartBlock(ModelPart& rMainModelPart, const std::string& InitialTabulation);

    void ReadSubModelPartDataBlock(ModelPart& rModelPart);

    void ReadSubModelPartTablesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart);

    void ReadSubModelPartPropertiesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart);

    /**
     * @brief Reads the nodes block of a submodelpart.
     * @param rMainModelPart Reference to the main model part.
     * @param rSubModelPart Reference to the submodelpart where the nodes will be added.
     */
    void ReadSubModelPartNodesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart);

    /**
     * @brief Reads the elements block of a submodelpart.
     * @param rMainModelPart Reference to the main model part.
     * @param rSubModelPart Reference to the submodelpart where the elements will be added.
     */
    void ReadSubModelPartElementsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart);

    /**
     * @brief Reads the conditions block of a submodelpart.
     * @param rMainModelPart Reference to the main model part.
     * @param rSubModelPart Reference to the submodelpart where the conditions will be added.
     */
    void ReadSubModelPartConditionsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart);

    /**
     * @brief Reads the geometries block of a submodelpart.
     * @param rMainModelPart Reference to the main model part.
     * @param rSubModelPart Reference to the submodelpart where the geometries will be added.
     */
    void ReadSubModelPartGeometriesBlock(
        ModelPart &rMainModelPart,
        ModelPart &rSubModelPart);

    /**
     * @brief Reads the master-slave constraints block of a submodelpart.
     * @param rMainModelPart Reference to the main model part.
     * @param rSubModelPart Reference to the submodelpart where the master-slave constraints will be added.
     */
    void ReadSubModelPartMasterSlaveConstraintsBlock(
        ModelPart& rMainModelPart,
        ModelPart& rSubModelPart
        );

    void DivideInputToPartitionsImpl(
        OutputFilesContainerType& rOutputFiles,
        SizeType NumberOfPartitions,
        const PartitioningInfo& rPartitioningInfo);

    void DivideModelPartDataBlock(OutputFilesContainerType& OutputFiles);

    void DivideTableBlock(OutputFilesContainerType& OutputFiles);

    void DividePropertiesBlock(OutputFilesContainerType& OutputFiles);

    /**
     * @brief Divides the nodes block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param NodesAllPartitions The partition indices for all nodes.
     */
    void DivideNodesBlock(OutputFilesContainerType& OutputFiles,
                          PartitionIndicesContainerType const& NodesAllPartitions);

    /**
     * @brief Divides the geometries block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param GeometriesAllPartitions The partition indices for all geometries.
     */
    void DivideGeometriesBlock(OutputFilesContainerType& OutputFiles,
                             PartitionIndicesContainerType const& GeometriesAllPartitions);

    /**
     * @brief Divides the elements block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ElementsAllPartitions The partition indices for all elements.
     */
    void DivideElementsBlock(OutputFilesContainerType& OutputFiles,
                             PartitionIndicesContainerType const& ElementsAllPartitions);

    /**
     * @brief Divides the conditions block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ConditionsAllPartitions The partition indices for all conditions.
     */
    void DivideConditionsBlock(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& ConditionsAllPartitions);

    /**
     * @brief Divides the master-slave constraints block into partitions.
     * @param rOutputFiles The container of output files for each partition.
     * @param rMasterSlaveConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
    void DivideMasterSlaveConstraintsBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

    void DivideNodalDataBlock(OutputFilesContainerType& OutputFiles,
                              PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideFlagVariableData(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideDofVariableData(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideScalarVariableData(OutputFilesContainerType& OutputFiles,
                                PartitionIndicesContainerType const& EntitiesPartitions,
                                std::string BlockName);

    template<class TValueType>
    void DivideVectorialVariableData(OutputFilesContainerType& OutputFiles,
                                     PartitionIndicesContainerType const& EntitiesPartitions,
                                     std::string BlockName);


    /**
     * @brief Divides the elemental data block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ElementsAllPartitions The partition indices for all elements.
     */
    void DivideElementalDataBlock(OutputFilesContainerType& OutputFiles,
                                  PartitionIndicesContainerType const& ElementsAllPartitions);

    /**
     * @brief Divides the conditional data block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ConditionsAllPartitions The partition indices for all conditions.
     */
    void DivideConditionalDataBlock(OutputFilesContainerType& OutputFiles,
                                    PartitionIndicesContainerType const& ConditionsAllPartitions);

    /**
     * @brief Divides the master-slave constraint data block into partitions.
     * @param rOutputFiles The container of output files for each partition.
     * @param rMasterSlaveConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
    void DivideMasterSlaveConstraintDataBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

    /**
     * @brief Divides the mesh block into partitions.
     * @param rOutputFiles The container of output files for each partition.
     * @param rNodesAllPartitions The partition indices for all nodes.
     * @param rElementsAllPartitions The partition indices for all elements.
     * @param rConditionsAllPartitions The partition indices for all conditions.
     * @param rMasterSlaveConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
    void DivideMeshBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rNodesAllPartitions,
        const PartitionIndicesContainerType& rElementsAllPartitions,
        const PartitionIndicesContainerType& rConditionsAllPartitions,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

    /**
     * @brief Divides the submodelpart block into partitions.
     * @param rOutputFiles The container of output files for each partition.
     * @param rNodesAllPartitions The partition indices for all nodes.
     * @param rElementsAllPartitions The partition indices for all elements.
     * @param rConditionsAllPartitions The partition indices for all conditions.
     * @param rMasterSlaveConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
	void DivideSubModelPartBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rNodesAllPartitions,
        const PartitionIndicesContainerType& rElementsAllPartitions,
        const PartitionIndicesContainerType& rConditionsAllPartitions,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

    void DivideMeshDataBlock(OutputFilesContainerType& OutputFiles);

    /**
     * @brief Divides the mesh nodes block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param NodesAllPartitions The partition indices for all nodes.
     */
    void DivideMeshNodesBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& NodesAllPartitions);

    /**
     * @brief Divides the mesh elements block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ElementsAllPartitions The partition indices for all elements.
     */
    void DivideMeshElementsBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& ElementsAllPartitions);

    /**
     * @brief Divides the mesh conditions block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ConditionsAllPartitions The partition indices for all conditions.
     */
    void DivideMeshConditionsBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions);

    /**
     * @brief Divides the mesh master-slave constraints block into partitions.
     * @param rOutputFiles The container of output files for each partition.
     * @param rMasterSlaveConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
    void DivideMeshMasterSlaveConstraintsBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

    /**
     * @brief Divides the submodelpart data block into partitions.
     * @param OutputFiles The container of output files for each partition.
     */
    void DivideSubModelPartDataBlock(OutputFilesContainerType& OutputFiles);

    /**
     * @brief Divides the submodelpart table block into partitions.
     * @param OutputFiles The container of output files for each partition.
     */
    void DivideSubModelPartTableBlock(OutputFilesContainerType& OutputFiles);

    /**
     * @brief Divides the submodelpart nodes block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param NodesAllPartitions The partition indices for all nodes.
     */
    void DivideSubModelPartNodesBlock(OutputFilesContainerType& OutputFiles,
        PartitionIndicesContainerType const& NodesAllPartitions);

    /**
     * @brief Divides the submodelpart elements block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ElementsAllPartitions The partition indices for all elements.
     */
    void DivideSubModelPartElementsBlock(OutputFilesContainerType& OutputFiles,
        PartitionIndicesContainerType const& ElementsAllPartitions);

    /**
     * @brief Divides the submodelpart conditions block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ConditionsAllPartitions The partition indices for all conditions.
     */
    void DivideSubModelPartConditionsBlock(OutputFilesContainerType& OutputFiles,
        PartitionIndicesContainerType const& ConditionsAllPartitions);

    /**
     * @brief Divides the submodelpart master-slave constraints block into partitions.
     * @param OutputFiles The container of output files for each partition.
     * @param ConstraintsAllPartitions The partition indices for all master-slave constraints.
     */
    void DivideSubModelPartMasterSlaveConstraintsBlock(
        OutputFilesContainerType& rOutputFiles,
        const PartitionIndicesContainerType& rMasterSlaveConstraintsAllPartitions
        );

	void WritePartitionIndices(OutputFilesContainerType& OutputFiles, PartitionIndicesType const&  NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions);

    void WriteCommunicatorData(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                               PartitionIndicesType const& NodesPartitions,
                               PartitionIndicesType const& ElementsPartitions,
                               PartitionIndicesType const& ConditionsPartitions,
                               PartitionIndicesContainerType const& NodesAllPartitions,
                               PartitionIndicesContainerType const& ElementsAllPartitions,
                               PartitionIndicesContainerType const& ConditionsAllPartitions);

    void WriteCommunicatorLocalNodes(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, PartitionIndicesType const& NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions);

    void WriteInAllFiles(OutputFilesContainerType& OutputFiles, std::string const& ThisWord);


    template<class TContainerType, class TKeyType>
    typename TContainerType::iterator FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName);

    // Basically it starts to read the character sequence until reaching a
    // "(" and then goes until corresponding ")" which means the vector or
    // matrix value is completely read. It can be used to read any kind of
    // vector or matrix with operator >> defined and writtern in following
    // format for a vector: [size] ( value1, value2,...., valueN )
    // format for a matrix: [size1,size2] ( )( )...( ) //look props read
    template<class TValueType>
    TValueType& ReadVectorialValue(TValueType& rValue);

    template<class TValueType>
    TValueType& ExtractValue(std::string rWord, TValueType & rValue);

    bool& ExtractValue(std::string rWord, bool & rValue);

    void ReadConstitutiveLawValue(ConstitutiveLaw::Pointer& rValue);

    ModelPartIO& ReadWord(std::string& Word);

    ModelPartIO& ReadBlock(std::string& Block, std::string const& BlockName);

    char SkipWhiteSpaces();

    char GetCharacter();

    void CheckStatement(std::string const& rStatement, std::string const& rGivenWord) const;

    void ResetInput();

    inline void CreatePartition(unsigned int NumberOfThreads,const int NumberOfRows, DenseVector<unsigned int>& partitions);

    /// Iterate over a Node block, calling ReorderedNodeId on each node.
    /** This method allows derived implementations to initialize reordering
     *  without storing the nodes.
     */
    void ScanNodeBlock();

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

}; // Class ModelPartIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  ModelPartIO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const ModelPartIO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.
