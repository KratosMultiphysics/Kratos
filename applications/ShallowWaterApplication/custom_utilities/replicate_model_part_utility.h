//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_REPLICATE_MODEL_PART_UTILITY
#define KRATOS_REPLICATE_MODEL_PART_UTILITY


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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
 * @ingroup ShallowWaterApplication
 * @class ReplicateModelPartUtility
 * @brief This utility replicates a model part to print the topography in the post-process
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) ReplicateModelPartUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Pointer definition of ReplicateModelPartUtility
    KRATOS_CLASS_POINTER_DEFINITION(ReplicateModelPartUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ReplicateModelPartUtility(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, bool ReplicateSubModelParts = false);

    /// Destructor.
    ~ReplicateModelPartUtility() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief It creates a copy of the origin model part and all the nodes, elements, conditions, properties, it also copies the ProcessInfo.
     */
    void Replicate();

    /**
     * @brief This method copies a variable from the origin model part to the destination model part for post-process purpose.
     */
    template<class TVarType>
    void TransferVariable(const TVarType& rVariable)
    {
        block_for_each(mrOriginModelPart.Nodes(), [&](NodeType& rNode){
            const auto dest_node = mReplicatedNodesMap[rNode.Id()];
            dest_node->FastGetSolutionStepValue(rVariable) = rNode.FastGetSolutionStepValue(rVariable);
        });
    }

    /**
     * @brief This method copies a variable from the origin model part to the destination model part for post-process purpose.
     */
    template<class TVarType>
    void TransferNonHistoricalVariable(const TVarType& rVariable)
    {
        block_for_each(mrOriginModelPart.Nodes(), [&](NodeType& rNode){
            const auto dest_node = mReplicatedNodesMap[rNode.Id()];
            dest_node->SetValue(rVariable, rNode.GetValue(rVariable));
        });
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
    std::string Info() const
    {
        return "ReplicateModelPartUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    const bool mReplicateSubModelParts;
    std::unordered_map<IndexType, Node<3>::Pointer> mReplicatedNodesMap;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void GetMaximumIds(IndexType& rUniqueNodeId, IndexType& rUniqueElemId, IndexType& rUniqueCondId, IndexType& rUniquePropId);

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
    ReplicateModelPartUtility& operator=(ReplicateModelPartUtility const& rOther) = delete;

    /// Copy constructor.
    ReplicateModelPartUtility(ReplicateModelPartUtility const& rOther) = delete;


    ///@}

}; // Class ReplicateModelPartUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ReplicateModelPartUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ReplicateModelPartUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REPLICATE_MODEL_PART_UTILITY  defined
