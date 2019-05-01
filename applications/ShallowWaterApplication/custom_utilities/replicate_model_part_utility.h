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
#include "includes/define.h"
#include "includes/model_part.h"


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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) ReplicateModelPartUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

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

    void Replicate();

    template<class TVarType>
    void TransferVariable(TVarType& rVariable)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfNodes()); ++i)
        {
            const auto it_node = mrOriginModelPart.NodesBegin() + i;
            const auto dest_node = mReplicatedNodes[it_node->Id()];
            dest_node->FastGetSolutionStepValue(rVariable) = it_node->FastGetSolutionStepValue(rVariable);
        }
    }

    template<class TVarType>
    void TransferNonHistoricalVariable(TVarType& rVariable)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfNodes()); ++i)
        {
            const auto it_node = mrOriginModelPart.NodesBegin() + i;
            const auto dest_node = mReplicatedNodes[it_node->Id()];
            dest_node->SetValue(rVariable, it_node->GetValue(rVariable));
        }
    }

    void SetOriginMeshPosition();

    void SetOriginMeshPosition(Variable<double>& rVariable);

    void SetDestinationMeshPosition();

    void SetDestinationMeshPosition(Variable<double>& rVariable);

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

    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    bool mReplicateSubModelParts;
    std::unordered_map<IndexType, Node<3>::Pointer> mReplicatedNodes;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void GetMaximumIds(IndexType& rUniqueNodeId, IndexType& rUniqueElemId, IndexType& rUniqueCondId);

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
    ReplicateModelPartUtility& operator=(ReplicateModelPartUtility const& rOther);

    /// Copy constructor.
    ReplicateModelPartUtility(ReplicateModelPartUtility const& rOther);


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
