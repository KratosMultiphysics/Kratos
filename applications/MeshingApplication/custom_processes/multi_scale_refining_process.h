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
#include "custom_utilities/uniform_refine_utility.h"


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
 * This process manages two model parts, the origin or the coarse model part
 * and the refined or the subscale model part
 * This process can be constructed again with the subscale as the origin model part
 * to get several subscales levels
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

    /**
     * Maps for AssignUniqueModelPartCollectionTagUtility
     */
    typedef std::unordered_map<IndexType, IndexType> IndexIndexMapType;
    typedef std::unordered_map<IndexType, std::vector<std::string>> IndexStringMapType;
    typedef std::unordered_map<IndexType, std::vector<IndexType>> IndexVectorMapType;

    /**
     * Map type to locate created refining interface conditions
     */
    typedef std::set<IndexVectorType> IndexVectorSetType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MultiScaleRefiningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MultiScaleRefiningProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    MultiScaleRefiningProcess(
        ModelPart& rThisCoarseModelPart,
        ModelPart& rThisRefinedModelPart,
        ModelPart& rThisVisualizationModelPart,
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

    /**
     * ExecuteRefinement transfers the entities TO_REFINE from the coarse to the refined model part and executes the refinement
     */
    void ExecuteRefinement();

    /**
     * ExecuteCoarsening deletes the entities of the refined model part which in the coarse model part are not TO_REFINE
     */
    void ExecuteCoarsening();

    /**
     * @brief InitializeNewModelPart is an auxiliary function to copy properties, variables, tables and sub model parts
     * @param rReferenceModelPart
     * @param rNewModelPart
     */
    static void InitializeNewModelPart(ModelPart& rReferenceModelPart, ModelPart& rNewModelPart);


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

    ModelPart& mrCoarseModelPart;  /// The coarse sub model part
    ModelPart& mrRefinedModelPart; /// Where the refinement is performed
    ModelPart& mrVisualizationModelPart;
    Parameters mParameters;

    unsigned int mEchoLevel;
    unsigned int mDivisionsAtSubscale;

    UniformRefineUtility<2> mUniformRefinement; /// The utility to perform the refinement

    IndexNodeMapType mCoarseToRefinedNodesMap; /// Mapping from coarse to refined
    IndexNodeMapType mRefinedToCoarseNodesMap; /// Mapping from refined to coarse

    std::string mRefinedInterfaceName;
    std::string mInterfaceConditionName;
    IndexVectorSetType mCoarseInterfacesSet;

    IndexStringMapType mCollections;  /// For AssignUniqueModelCollectionTagUtility

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Perform a check with the parameters
     */
    void Check();

    /**
     * @brief
     */
    void InterpolateLevelBoundaryValuesAtSubStep(const int& rSubStep, const int& rSubSteps);

    /**
     * @brief
     */
    void UpdateSubLevel();

    /**
     * @brief
     */
    void TransferDataToCoarseLevel();

    /* TODO: debug and move to AssignUniqueModelPArtCollectionTagUtility */
    StringVectorType RecursiveGetSubModelPartNames(ModelPart& rThisModelPart, std::string Prefix = "");

    /* TODO: debug and move to AssignUniqueModelPArtCollectionTagUtility */
    ModelPart& RecursiveGetSubModelPart(ModelPart& rThisModelPart, std::string FullName);

    /**
     * @brief InitializeCoarseModelPart
     * @param rNames Is the vector containing the sub model part names
     */
    void InitializeCoarseModelPart();

    /**
     * @brief InitializeRefinedModelPart creates the refined sub model part
     * @detail The method copy the model part hierarchy from the coarse to the refined model part
     * @param rNames The vector containing the sub model part names
     */
    void InitializeRefinedModelPart(const StringVectorType& rNames);

    /**
     * @brief InitializeVisualizationModelPart adds all the nodes, elements and conditions to the visualization model part
     * @param rNemes The vector containing the sub model part names
     */
    void InitializeVisualizationModelPart(const StringVectorType& rNames); 

    /**
     * @brief AddVariablesToRefinedModelPart copies the variables from the coarse to the refined model part
     */
    void AddVariablesToRefinedModelPart();

    /**
     * @brief AddAllPropertiesToModelPart adds all properties from an origin model part to a destination model part
     * @param rOriginModelPart
     * @param pDestinationModelPart
     */
    static void AddAllPropertiesToModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart);

    /**
     * @brief AddAllTablesToModelPart adds all tables from an origin model part to a destination model part
     * @param rOriginModelPart
     * @param pDestinationModelPart
     */
    static void AddAllTablesToModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart);

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
    void CreateElementsToRefine(IndexType& rElemId, IndexIndexMapType& rElemTag);

    /**
     * @brief Create the auxiliary conditions in the refined sub model part
     * @param rCondId the condition Id will be ++rCondId
     */
    void CreateConditionsToRefine(IndexType& rCondId, IndexIndexMapType& rCondTag);

    /**
     * This method makes a substitution on the visualization model part: removes the coarse entities and adds the refined ones
     */
    void UpdateVisualizationAfterRefinement();

    /**
     * @brief IdentifyNodesToErase looks for the nodes which should be removed from the refined model part
     * @detail When a node is not TO_REFINE and is currently refined, sets OLD_ENTITY flag in the coarse
     * model part and remove it from the unordered_maps
     * @see CloneNodesToRefine
     */
    void IdentifyParentNodesToErase();

    /**
     * @brief IdentifyElementsToErase looks for the elements which should be removed from the refined model part
     * @detail Sets TO_ERASE flag in the refined model part when a node in the coarse model part is OLD_ENTITY
     * @see IdentifyParentNodesToErase
     */
    void IdentifyElementsToErase();

    /**
     * @brief IdentifyConditionsToErase looks for the condtions which should be removed from the refined model part
     * @detail Sets TO_ERASE flag in the refined model part when a node in the coarse model part is OLD_ENTITY
     * @see IdentifyParentNodesToRefine
     */
    void IdentifyConditionsToErase();

    /**
     * @brief IdentifyrefinedNodesToErase looks for the nodes which should be removed from the refiend model part
     * @detail When a refined element is TO_ERASE, sets TO_ERASE flag to all its nodes.
     * If a node is TO_REFINE (e.g. the refining interface), the element is not TO_ERASE
     * @see IdentifyElementsToErase
     */
    void IdentifyRefinedNodesToErase();

    /**
     * @brief FinalizeRefinement resets the flags on the nodes and elements
     * @detail NEW_ENTITY is set to false
     * @see CloneNodesToRefine
     */
    void FinalizeRefinement();

    /**
     * @brief FinalizeCoarsening resets the flags on the nodes, elements and conditions
     * @detail MeshingFlags::TO_COARSE
     */
    void FinalizeCoarsening();

    /**
     * @brief IdentifyRefiningInterface
     * TODO: define how it works
     */
    void IdentifyRefiningInterface(IndexType& rCondId);

    /**
     * @brief GetLastId gets the absolute root model part and looks for the maximum id's
     */
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
