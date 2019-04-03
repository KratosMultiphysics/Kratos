//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                    
//

#if !defined(KRATOS_FAST_TRANSFER_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED )
#define  KRATOS_FAST_TRANSFER_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"

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
        
///@name Kratos Classes
///@{

/**
 * @class FastTransferBetweenModelPartsProcess
 * @ingroup KratosCore
 * @brief The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
 * @details This function assigns a value to a variable belonging to all of the nodes in a given mesh
 * @note If the MODIFIED flag is set it will create new entities, copying the geometry and data, instead of just transfering
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) FastTransferBetweenModelPartsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The type used for the node
    typedef Node<3> NodeType;

    // General containers type definitions
    typedef ModelPart::NodesContainerType                                 NodesArrayType;
    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;
    typedef ModelPart::ElementsContainerType                           ElementsArrayType;
    typedef ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintArrayType;

    // General containers iterators type definitions
    typedef NodesArrayType::iterator                                  IteratorNodesArrayType;
    typedef ConditionsArrayType::iterator                        IteratorConditionsArrayType;
    typedef ElementsArrayType::iterator                            IteratorElementsArrayType;
    typedef MasterSlaveConstraintArrayType::iterator IteratorMasterSlaveConstraintsArrayType;

    /// The type used to identify the size
    typedef std::size_t SizeType;

    /// Pointer definition of FastTransferBetweenModelPartsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FastTransferBetweenModelPartsProcess);

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enum helps us to identify the elements to transfer between the modelparts
     * @todo Add intermediate combinations of constraints, conditions and elements
     */
    enum class EntityTransfered {
        NODES = 0,
        ELEMENTS = 1,
        NODESANDELEMENTS = 2,
        CONDITIONS = 3,
        NODESANDCONDITIONS = 4,
        CONSTRAINTS = 5,
        NODESANDCONSTRAINTS = 6,
        ALL = 7
    };

    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Default constructor. Without flag
     * @param rDestinationModelPart The destination model part
     * @param rOriginModelPart The origin model part
     * @param Entity The elements to transfer
     * @param Flag The flag used to differentiate between elements to transfer
     * @param ReplicateEntities  If the entities are replicated or transfered
     */
    FastTransferBetweenModelPartsProcess(
        ModelPart& rDestinationModelPart,
        ModelPart& rOriginModelPart,
        const EntityTransfered Entity = EntityTransfered::ALL,
        const Flags Flag = Flags(),
        const bool ReplicateEntities = false
        );

    /// Destructor.
    ~FastTransferBetweenModelPartsProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()();


    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the FastTransferBetweenModelPartsProcess algorithms.
    void Execute() override;

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
        return "FastTransferBetweenModelPartsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FastTransferBetweenModelPartsProcess";
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

//     /// Copy constructor.
//     FastTransferBetweenModelPartsProcess(FastTransferBetweenModelPartsProcess const& rOther);

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

    ModelPart& mrDestinationModelPart; /// The destination model part
    
    ModelPart& mrOriginModelPart;      /// The origin model part
    
    const EntityTransfered mEntity;    /// The entity to transfer

    const Flags mFlag;                 /// A flag in order to tranfer only components with that flag
    
    ///@}
    ///@name Private Operators
    ///@{
    
    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief This method transfer the entities without considering the flags
     */
    void TransferWithoutFlags();

    /**
     * @brief This method transfer the entities considering the flags
     */
    void TransferWithFlags();

    /**
     * @brief This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     */

    void ReorderAllIds(ModelPart& rThhisModelPart);

    /**
     * @brief This method replicates the entities without considering the flags
     */
    void ReplicateWithoutFlags();

    /**
     * @brief This method replicates the entities considering the flags
     */
    void ReplicateWithFlags();

    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    FastTransferBetweenModelPartsProcess& operator=(FastTransferBetweenModelPartsProcess const& rOther);

    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class FastTransferBetweenModelPartsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FastTransferBetweenModelPartsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FastTransferBetweenModelPartsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FAST_TRANSFER_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED  defined
