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

/// IO provides different implementation of input output procedures which can be used to read and write with different formats and characteristics.
/** IO provides different implementation of input output procedures which can be used to read and write with different formats and characteristics.
 * An automatic configurable IO module is added to these components providing the complete set of solutions necessary for dealing with multi-disciplinary problems.
 * This IO module uses different component lists to adjust itself when reading and writing new concepts originating from different fields of analysis.
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

    virtual bool ReadNode(NodeType& rThisNode)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class." << std::endl;
    }

    virtual bool ReadNodes(NodesContainerType& rThisNodes)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual std::size_t ReadNodesNumber()
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class." << std::endl;;
    }

    virtual void WriteNodes(NodesContainerType const& rThisNodes)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadProperties(Properties& rThisProperties)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadProperties(PropertiesContainerType& rThisProperties)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteProperties(Properties const& rThisProperties)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteProperties(PropertiesContainerType const& rThisProperties)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual std::size_t  ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteElements(ElementsContainerType const& rThisElements)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual std::size_t  ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteConditions(ConditionsContainerType const& rThisConditions)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadInitialValues(ModelPart& rThisModelPart)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadInitialValues(NodesContainerType& rThisNodes, ElementsContainerType& rThisElements, ConditionsContainerType& rThisConditions)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

    virtual void ReadMesh(MeshType & rThisMesh)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void ReadModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void WriteNodeMesh( MeshType& rThisMesh )
    {
        KRATOS_ERROR <<  "Calling base class method (WriteNodeMesh). Please check the implementation of derived classes" << std::endl;
    }

    virtual std::size_t ReadNodalGraph(ConnectivitiesContainerType& aux_connectivities)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;;
    }

    virtual void DivideInputToPartitions(SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                                         PartitionIndicesType const& NodesPartitions,
                                         PartitionIndicesType const& ElementsPartitions,
                                         PartitionIndicesType const& ConditionsPartitions,
                                         PartitionIndicesContainerType const& NodesAllPartitions,
                                         PartitionIndicesContainerType const& ElementsAllPartitions,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
    }

    virtual void DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * Streams,
                                         SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                                         PartitionIndicesType const& NodesPartitions,
                                         PartitionIndicesType const& ElementsPartitions,
                                         PartitionIndicesType const& ConditionsPartitions,
                                         PartitionIndicesContainerType const& NodesAllPartitions,
                                         PartitionIndicesContainerType const& ElementsAllPartitions,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions)
    {
        KRATOS_ERROR <<  "Calling base class member. Please check the definition of derived class" << std::endl;
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
