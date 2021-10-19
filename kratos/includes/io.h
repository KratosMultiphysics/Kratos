//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_IO_H_INCLUDED )
#define  KRATOS_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <unordered_set>

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

    typedef DenseMatrix<int> GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IO() {}

    /// Destructor.
    virtual ~IO() {}


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

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

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
    virtual void WriteMesh( MeshType& rThisMesh )
    {
        KRATOS_ERROR << "Calling base class method (WriteMesh). Please check the implementation of derived classes" << std::endl;
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
    virtual void WriteModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_ERROR << "Calling base class method (WriteModelPart). Please check the definition of derived class" << std::endl;
    }

    /**
     * @brief This method writes the node mesh
     * @param rThisMesh The mesh to be written
     */
    virtual void WriteNodeMesh( MeshType& rThisMesh )
    {
        KRATOS_ERROR << "Calling base class method (WriteNodeMesh). Please check the implementation of derived classes" << std::endl;
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
     * @param rDomainsColoredGraph The colors of the partition graph
     * @param rNodesPartitions The partitions indices of the nodes
     * @param rElementsPartitions The partitions indices of the elements
     * @param rConditionsPartitions The partitions indices of the conditions
     * @param rNodesAllPartitions The partitions of the nodes
     * @param rElementsAllPartitions The partitions of the elements
     * @param rConditionsAllPartitions The partitions of the conditions
     */
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
     * @param rDomainsColoredGraph The colors of the partition graph
     * @param rNodesPartitions The partitions indices of the nodes
     * @param rElementsPartitions The partitions indices of the elements
     * @param rConditionsPartitions The partitions indices of the conditions
     * @param rNodesAllPartitions The partitions of the nodes
     * @param rElementsAllPartitions The partitions of the elements
     * @param rConditionsAllPartitions The partitions of the conditions
     */
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

    virtual void ReadSubModelPartElementsAndConditionsIds(
        std::string const& rModelPartName,
        std::unordered_set<SizeType> &rElementsIds,
        std::unordered_set<SizeType> &rConditionsIds)
    {
        KRATOS_ERROR << "Calling base class method (ReadSubModelPartElementsAndConditionsIds). Please check the definition of derived class" << std::endl;
    }

    virtual std::size_t ReadNodalGraphFromEntitiesList(
        ConnectivitiesContainerType& rAuxConnectivities,
        std::unordered_set<SizeType> &rElementsIds,
        std::unordered_set<SizeType> &rConditionsIds)
    {
        KRATOS_ERROR << "Calling base class method (ReadNodalGraphFromEntitiesList). Please check the definition of derived class" << std::endl;
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

#endif // KRATOS_IO_H_INCLUDED  defined
