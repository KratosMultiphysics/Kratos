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

    #if !defined(ENTITY_TRANSFERED)
    #define ENTITY_TRANSFERED
        enum EntityString {Nodes = 0, Elements = 1, NodesAndElements = 2, Conditions = 3, All = 4};
    #endif

///@}
///@name  Functions
///@{
        
///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class FastTransferBetweenModelPartsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

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
        
        const int num_nodes = mrOriginModelPart.Nodes().size();

        if (num_nodes != 0 && (mEntity == All || mEntity == Nodes || mEntity == NodesAndElements))
            mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

        const int num_elements = mrOriginModelPart.Elements().size();

        if (num_elements != 0 && (mEntity == All || mEntity == Elements || mEntity == NodesAndElements))
            mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

        const int num_conditions = mrOriginModelPart.Conditions().size();

        if (num_conditions != 0 && (mEntity == All || mEntity == Conditions))
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

    ModelPart& mrDestinationModelPart; // The destination model part
    
    ModelPart& mrOriginModelPart;      // The origin model part
    
    const EntityString mEntity;        // The entity to transfer
    
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

    EntityString ConvertEntity(const std::string& str)
    {
        if(str == "Nodes")
            return Nodes;
        else if(str == "Elements")
            return Elements;
        else if(str == "NodesAndElements")
            return NodesAndElements;
        else if(str == "Conditions")
            return Conditions;
        else if(str == "All")
            return All;
        else
            return All;
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
