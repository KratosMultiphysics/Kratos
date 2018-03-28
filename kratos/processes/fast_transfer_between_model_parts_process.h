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
#include "includes/kratos_parameters.h"
#include "processes/process.h"

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

    /**
     * @brief This enum helps us to identify the elements to transfer between the modelparts
     */
    enum class EntityTransfered {
        NODES = 0,
        ELEMENTS = 1,
        NODESANDELEMENTS = 2,
        CONDITIONS = 3,
        ALL = 4};

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
 * @author Vicente Mataix Ferrandiz
*/
class FastTransferBetweenModelPartsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The type used to identify the size
    typedef std::size_t SizeType;

    /// Pointer definition of FastTransferBetweenModelPartsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FastTransferBetweenModelPartsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    
    FastTransferBetweenModelPartsProcess(
        ModelPart& rDestinationModelPart,
        ModelPart& rOriginModelPart,
        const std::string EntityString
        ) : Process(),
            mrDestinationModelPart(rDestinationModelPart), 
            mrOriginModelPart(rOriginModelPart),
            mEntity(ConvertEntity(EntityString))
    {
        KRATOS_TRY
                
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~FastTransferBetweenModelPartsProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the FastTransferBetweenModelPartsProcess algorithms.
    void Execute() override
    {
        KRATOS_TRY;
        
        const SizeType num_nodes = mrOriginModelPart.Nodes().size();

        if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS))
            mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

        const SizeType num_elements = mrOriginModelPart.Elements().size();

        if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
            mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

        const SizeType num_conditions = mrOriginModelPart.Conditions().size();

        if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS))
             mrDestinationModelPart.AddConditions(mrOriginModelPart.ConditionsBegin(),mrOriginModelPart.ConditionsEnd());

        KRATOS_CATCH("");
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
    
    ///@}
    ///@name Private Operators
    ///@{
    
    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This converts the entity string to an enum
     * @param str The string that you want to convert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */

    EntityTransfered ConvertEntity(const std::string& str)
    {
        if(str == "Nodes")
            return EntityTransfered::NODES;
        else if(str == "Elements")
            return EntityTransfered::ELEMENTS;
        else if(str == "NodesAndElements")
            return EntityTransfered::NODESANDELEMENTS;
        else if(str == "Conditions")
            return EntityTransfered::CONDITIONS;
        else if(str == "All")
            return EntityTransfered::ALL;
        else
            return EntityTransfered::ALL;
    }
    
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
