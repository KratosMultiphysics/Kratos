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

    /// The type used for the node
    typedef Node<3> NodeType;

    /// The type used to identify the size
    typedef std::size_t SizeType;

    /// Pointer definition of FastTransferBetweenModelPartsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FastTransferBetweenModelPartsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Default constructor. Without flag
     * @param rDestinationModelPart The destination model part
     * @param rOriginModelPart The origin model part
     * @param EntityString The elements to transfer
     */
    FastTransferBetweenModelPartsProcess(
        ModelPart& rDestinationModelPart,
        ModelPart& rOriginModelPart,
        const std::string EntityString = "All",
        const std::string FlagName = ""
        ) : Process(),
            mrDestinationModelPart(rDestinationModelPart), 
            mrOriginModelPart(rOriginModelPart),
            mEntity(ConvertEntity(EntityString)),
            mFlagName(FlagName)
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
        
        // In case of not flag defined we transfer all the elements
        if (!KratosComponents<Flags>::Has(mFlagName)) {
            const SizeType num_nodes = mrOriginModelPart.Nodes().size();

            if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS))
                mrDestinationModelPart.AddNodes(mrOriginModelPart.NodesBegin(),mrOriginModelPart.NodesEnd());

            const SizeType num_elements = mrOriginModelPart.Elements().size();

            if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
                mrDestinationModelPart.AddElements(mrOriginModelPart.ElementsBegin(),mrOriginModelPart.ElementsEnd());

            const SizeType num_conditions = mrOriginModelPart.Conditions().size();

            if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS))
                mrDestinationModelPart.AddConditions(mrOriginModelPart.ConditionsBegin(),mrOriginModelPart.ConditionsEnd());
        } else {
//             const Flags this_flag = KratosComponents<Flags>::Get(mFlagName);
//
//             std::vector<NodeType::Pointer> vector_nodes;
//             std::vector<Element::Pointer> vector_elements;
//             std::vector<Condition::Pointer> vector_conditions;
//
//             // Creating a buffer for parallel vector fill
//             const int num_threads = OpenMPUtils::GetNumThreads();
//             std::vector<std::vector<NodeType::Pointer>> nodes_buffer(num_threads);
//             std::vector<std::vector<Element::Pointer>> elements_buffer(num_threads);
//             std::vector<std::vector<Condition::Pointer>> conditions_buffer(num_threads);
//
//             // Auxiliar sizes
//             const int num_nodes = static_cast<int>(mrOriginModelPart.Nodes().size());
//             const int num_elements = static_cast<int>(mrOriginModelPart.Elements().size());
//             const int num_conditions = static_cast<int>(mrOriginModelPart.Conditions().size());
//
//             #pragma omp parallel
//             {
//                 const int id = OpenMPUtils::ThisThread();
//
//                 if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS)) {
//                     #pragma omp for
//                     for(int i = 0; i < num_nodes; ++i) {
//                         auto it_node = mrOriginModelPart.NodesBegin() + i;
//                         if (it_node->Is(this_flag)) {
//                             (nodes_buffer[id]).push_back(*(it_node.base()));
//                         }
//                     }
//                 }
//
//                 if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS)) {
//                     #pragma omp for
//                     for(int i = 0; i < num_elements; ++i) {
//                         auto it_elem = mrOriginModelPart.ElementsBegin() + i;
//                         if (it_elem->Is(this_flag)) {
//                             (elements_buffer[id]).push_back(*(it_elem.base()));
//                         }
//                     }
//                 }
//
//                 if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS)) {
//                     #pragma omp for
//                     for(int i = 0; i < num_nodes; ++i) {
//                         auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
//                         if (it_cond->Is(this_flag)) {
//                             (conditions_buffer[id]).push_back(*(it_cond.base()));
//                         }
//                     }
//                 }
//
//                 // Combine buffers together
//                 #pragma omp single
//                 {
//                     for( auto& node_buffer : nodes_buffer)
//                         std::move(node_buffer.begin(),node_buffer.end(),back_inserter(vector_nodes));
//                     for( auto& element_buffer : elements_buffer)
//                         std::move(element_buffer.begin(),element_buffer.end(),back_inserter(vector_elements));
//                     for( auto& condition_buffer : conditions_buffer)
//                         std::move(condition_buffer.begin(),condition_buffer.end(),back_inserter(vector_conditions));
//                 }
//             }
//
//             // We transfer
//             if (num_nodes != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::NODES || mEntity == EntityTransfered::NODESANDELEMENTS))
//                 mrDestinationModelPart.AddNodes(vector_nodes.begin(),vector_nodes.end());
//
//             if (num_elements != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::ELEMENTS || mEntity == EntityTransfered::NODESANDELEMENTS))
//                 mrDestinationModelPart.AddElements(vector_elements.begin(),vector_elements.end());
//
//             if (num_conditions != 0 && (mEntity == EntityTransfered::ALL || mEntity == EntityTransfered::CONDITIONS))
//                 mrDestinationModelPart.AddConditions(vector_conditions.begin(),vector_conditions.end());
        }

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

    const std::string mFlagName;       /// A flag in order to tranfer only components with that flag
    
    ///@}
    ///@name Private Operators
    ///@{
    
    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief This converts the entity string to an enum
     * @param Str The string that you want to convert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */

    EntityTransfered ConvertEntity(const std::string& Str)
    {
        if(Str == "Nodes" || Str == "NODES")
            return EntityTransfered::NODES;
        else if(Str == "Elements" || Str == "ELEMENTS")
            return EntityTransfered::ELEMENTS;
        else if(Str == "NodesAndElements" || Str == "NODESANDELEMENTS")
            return EntityTransfered::NODESANDELEMENTS;
        else if(Str == "Conditions" || Str == "CONDITIONS")
            return EntityTransfered::CONDITIONS;
        else if(Str == "All" || Str == "ALL")
            return EntityTransfered::ALL;
        else
            KRATOS_ERROR << "The entity declared " << Str << " is not on the following list\n" \
                         << "- NODES\n" << "- ELEMENTS\n" << "- NODESANDELEMENTS\n" << "- CONDITIONS\n" << "- ALL\n" << std::endl;
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
