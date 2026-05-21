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

#if !defined(KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED )
#define  KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part_io.h"


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
/** This class writes all modelpart data including the meshes.
*/
class KRATOS_API(KRATOS_CORE) ReorderConsecutiveModelPartIO : public ModelPartIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReorderConsecutiveModelPartIO
    KRATOS_CLASS_POINTER_DEFINITION(ReorderConsecutiveModelPartIO);

    typedef ModelPartIO BaseType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::NodesContainerType NodesContainerType;

    typedef BaseType::PropertiesContainerType PropertiesContainerType;

    typedef BaseType::ElementsContainerType ElementsContainerType;

    typedef BaseType::ConditionsContainerType ConditionsContainerType;

    typedef BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

    typedef BaseType::OutputFilesContainerType OutputFilesContainerType;

    typedef std::map<SizeType, SizeType> IdMapType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    // inherit constructors of baseclass
    using ModelPartIO::ModelPartIO;

    // this (inherited) constructor cannot be used with this class
    ReorderConsecutiveModelPartIO(
        Kratos::shared_ptr<std::iostream> Stream,
        const Flags Options) = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    SizeType mNumberOfNodes = 0;
    SizeType mNumberOfElements = 0;
    SizeType mNumberOfConditions = 0;
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    SizeType ReorderedNodeId(ModelPartIO::SizeType NodeId) override;
    SizeType ReorderedElementId(ModelPartIO::SizeType ElementId) override;
    SizeType ReorderedConditionId(ModelPartIO::SizeType ConditionId) override;

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

    ///@name Member Variables
    ///@{

    IdMapType mNodeIdMap = {};
    IdMapType mElementIdMap = {};
    IdMapType mConditionIdMap = {};

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

}; // Class ReorderConsecutiveModelPartIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  ReorderConsecutiveModelPartIO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const ReorderConsecutiveModelPartIO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.

#endif // KRATOS_REORDER_CONSECUTIVE_MODEL_PART_IO_H_INCLUDED  defined
