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

#if !defined( KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED )
#define KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos {

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

/// This class provides a refining utility to perform multi scale analysis
/**
 * The process creates a nested model part with the origin model part and
 * another nested model part with the refined mesh
 * This process can be called within the refined sub model part
 * The origin model part is stored under a new sub model part (own):
 * MainModelPart
 *     Own
 *     RefinedSubModelPart
 *         Own
 *         RefinedSubModelPart
 *             Own
 *             ...
 */
class MultiScaleRefiningProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /**
     * Node type
     */
    typedef Node<3> NodeType;

    /**
     * Type of IDs
     */
    typedef std::size_t IndexType;

    /**
     * Vector of IndexType
     */
    typedef std::vector<IndexType> IndexVectorType;

    /**
     * Vector of strings type
     */
    typedef std::vector<std::string> StringVectorType;

    /**
     * Map types to locate nodes in the mesh
     */
    typedef std::unordered_map<IndexType, NodeType::Pointer> IndexNodeMapType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MultiScaleRefiningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MultiScaleRefiningProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    MultiScaleRefiningProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~MultiScaleRefiningProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override {}

    void ExecuteRefinement();

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
        return "MultiScaleRefiningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    ModelPart& mrRootModelPart;
    Parameters mParameters;

    std::string mOwnName;     /// The coarse sub model part
    std::string mRefinedName; /// Where the refinement is performed
    std::string mElementName;
    std::string mConditionName;

    unsigned int mEchoLevel;
    unsigned int mDivisions;

    ModelPart::Pointer mpOwnModelPart;     /// The coarse sub model part
    ModelPart::Pointer mpRefinedModelPart; /// Where the refinement is performed

    StringVectorType mModelPartsNames;  /// The names of the sub model parts hierarchy

    IndexNodeMapType mCoarseToRefinedNodesMap; /// Mapping from own to refined
    IndexNodeMapType mRefinedToCoarseNodesMap; /// Mapping from refined to own

    std::string mInterfaceName;
    std::string mInterfaceConditionName;

    ///@}
    ///@name Private Operators
    ///@{

    void InterpolateLevelBoundaryValuesAtSubStep(const int& rSubStep, const int& rSubSteps);

    void UpdateSubLevel();

    void TransferDataToCoarseLevel();

    StringVectorType RecursiveGetSubModelPartNames(ModelPart& rThisModelPart, std::string Prefix = "");

    ModelPart& RecursiveGetSubModelPart(ModelPart& rThisModelPart, std::string FullName);

    void InitializeOwnModelPart(const StringVectorType& rNames);

    void InitializeOwnModelPart(const std::string& rOwnName, const StringVectorType& rNames);  // TODO: remove this method

    void InitializeRefinedModelPart(const StringVectorType& rNames);

    void InitializeRefinedModelPart(const std::string& rRefinedName, const std::string& rOwnName, const StringVectorType& rNames); // TODO: remove this method

    void AddAllPropertiesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart);

    void AddAllTablesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart);

    void AddAllNodesToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart);

    void AddAllElementsToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart);

    void AddAllConditionsToModelPart(ModelPart& rOriginModelPart, ModelPart::Pointer pDestinationModelPart);

    /**
     * @brief This function sets the elements TO_REFINE depending on the nodal flags
     * @detail An element is TO_REFINE if all the nodes are TO_REFINE and, at least one node is NEW_ENTITY
     * @see CloneNodesToRefine
     */
    void MarkElementsFromNodalFlag();

    /**
     * @brief This function sets the conditions TO_REFINE depending on the nodal flags
     * @detail An condition is TO_REFINE if all the nodes are TO_REFINE and, at least one node is NEW_ENTITY
     * @see CloneNodesToRefine
     */
    void MarkConditionsFromNodalFlag();

    /**
     * @brief This function creates a copy of the nodes on the refined sub model part
     * @detail Only are copied (NEW_ENTITY) the nodes which are not already present in the refined sub model part
     * @param rNodeId the node Id will be ++rNodeId
     */
    void CloneNodesToRefine(IndexType& rNodeId);

    /**
     * @brief Create the auxiliary nodes in the refined sub model part
     * @param rElemId the element Id will be ++rElemId
     */
    void CreateElementsToRefine(IndexType& rElemId);

    /**
     * @brief Create the auxiliary conditions in the refined sub model part
     * @param rCondId the condition Id will be ++rCondId
     */
    void CreateConditionsToRefine(IndexType& rCondId);

    void IdentifyRefiningInterface();

    void GetLastId(IndexType& rNodesId, IndexType& rElemsId, IndexType& rCondsId);

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
    // MultiScaleRefiningProcess& operator=(MultiScaleRefiningProcess const& rOther);

    /// Copy constructor.
    // MultiScaleRefiningProcess(MultiScaleRefiningProcess const& rOther);

    ///@}

}; // Class MultiScaleRefiningProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MultiScaleRefiningProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MultiScaleRefiningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED
