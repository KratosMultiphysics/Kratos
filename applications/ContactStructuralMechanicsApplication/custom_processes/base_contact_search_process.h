// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_processes/normal_gap_process.h"
#include "spatial_containers/point_object.h"
#include "custom_conditions/paired_condition.h"
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    using SizeType = std::size_t;

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
 * @class BaseContactSearchProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process has as objective to create the contact conditions.
 * @details The conditions that can be created are Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) BaseContactSearchProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Containers definition
    using NodesArrayType = ModelPart::NodesContainerType;
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Geometry type definition
    using GeometryType = Geometry<Node>; 

    /// Definition of the indextype
    using IndexType = std::size_t;

    /// Type definitions for the tree
    using PointType = PointObject<Condition>; /// Point type definition
    using PointTypePointer = typename PointType::Pointer; /// Point type pointer definition
    using PointVector = std::vector<PointTypePointer>; /// Point vector definition

    /// KDtree definitions
    using BucketType = Bucket<3, PointType, PointVector>; /// Bucket type definition
    using KDTree = Tree<KDTreePartition<BucketType>>; /// KDTree type definition

    /// The type of mapper considered
    using NormalGapProcessType = NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>; /// Normal gap process type definition

    /// The definition of zero tolerance
    static constexpr double GapThreshold = 2.0e-3;

    /// The definition of zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of BaseContactSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( BaseContactSearchProcess );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( INVERTED_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG( CREATE_AUXILIAR_CONDITIONS );
    KRATOS_DEFINE_LOCAL_FLAG( MULTIPLE_SEARCHS );
    KRATOS_DEFINE_LOCAL_FLAG( PREDEFINE_MASTER_SLAVE );
    KRATOS_DEFINE_LOCAL_FLAG( PURE_SLIP );

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief An enumeration of the different types of search trees that can be used in the contact search process.
     */
    enum class SearchTreeType {
        KdtreeInRadius = 0, ///< A kd-tree search within a given radius.
        KdtreeInBox = 1, ///< A kd-tree search within a given box.
        KdtreeInRadiusWithOBB = 2, ///< A kd-tree search within a given radius, using oriented bounding boxes (OBBs).
        KdtreeInBoxWithOBB = 3, ///< A kd-tree search within a given box, using OBBs.
        OctreeWithOBB = 4, ///< An octree search using OBBs.
        Kdop = 5 ///< A kd-tree search using discrete oriented polytopes (DOPs).
    };

    /**
     * @brief An enumeration of the possible results of a check during the contact search process.
     */
    enum class CheckResult {
        Fail = 0, ///< The check failed.
        AlreadyInTheMap = 1, ///< The check succeeded, and the element is already in the map.
        OK = 2 ///< The check succeeded, and the element was added to the map.
    };

    /**
     * @brief An enumeration of the different types of gap checks that can be performed during the contact search process.
     */
    enum class CheckGap {
        NoCheck = 0, ///< No gap check is performed.
        DirectCheck = 1, ///< A direct gap check is performed.
        MappingCheck = 2 ///< A mapping gap check is performed.
    };

    /**
     * @brief An enumeration of the different types of solutions that can be used in the contact search process.
     */
    enum class TypeSolution {
        NormalContactStress = 0, ///< A normal contact stress solution.
        ScalarLagrangeMultiplier = 1, ///< A scalar Lagrange multiplier solution.
        VectorLagrangeMultiplier = 2, ///< A vector Lagrange multiplier solution.
        FrictionlessPenaltyMethod = 3, ///< A frictionless penalty method solution.
        FrictionalPenaltyMethod = 4, ///< A frictional penalty method solution.
        OtherFrictionless = 5, ///< Another frictionless solution.
        OtherFrictional = 6 ///< Another frictional solution.
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part to be considered
     * @param ThisParameters The configuration parameters, it includes:
     *                       - The allocation considered in the search
     *                       - The factor considered to check if active or not
     *                       - The integration order considered
     *                       - The size of the bucket
     *                       - The proportion increased of the Radius/Bounding-box volume for the search
     *                       - TypeSearch: 0 means search in radius, 1 means search in box
     * @param pPairedProperties Properties of the pair
     * @todo Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @note Use an InterfacePreprocess object to create such a model part from a regular one:
     *          -# InterfaceMapper = InterfacePreprocess()
     *          -# InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    BaseContactSearchProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters =  Parameters(R"({})"),
        Properties::Pointer pPairedProperties = nullptr
        );

    /// Destructor.
    ~BaseContactSearchProcess() override = default;

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

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function initializes the ALM frictionless mortar conditions already created
     */
    void InitializeMortarConditions();

    /**
     * @brief This function clears the mortar conditions already created
     */
    virtual void ClearMortarConditions();

    /**
     * @brief This method checks that the contact model part is unique (so the model parts contain unique contact pairs)
     */
    virtual void CheckContactModelParts();

    /**
     * @brief This function creates a lists  points ready for the Mortar method
     */
    void CreatePointListMortar();

    /**
     * @brief This function updates a lists  points ready for the Mortar method
     */
    void UpdatePointListMortar();

    /**
     * @brief This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     */
    void UpdateMortarConditions();

    /**
     * @brief It checks the current mortar conditions
     */
    void CheckMortarConditions();

    /**
     * @brief It sets if the search is inverted
     */
    void InvertSearch();

    /**
     * @brief This resets the contact operators
     */
     virtual void ResetContactOperators();

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns a string with information about the contact search process.
     * @return A string with information about the contact search process.
     */
    std::string Info() const override
    {
        return "BaseContactSearchProcess";
    }

    /**
     * @brief Prints information about the contact search process to an output stream.
     * @param rOStream The output stream to print to.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    ModelPart& mrMainModelPart;                            /// The main model part
    Parameters mThisParameters;                            /// The configuration parameters
    CheckGap mCheckGap;                                    /// If the gap is checked during the search
    TypeSolution mTypeSolution;                            /// The solution type
    PairedCondition const* mpReferenceCondition = nullptr; /// The prototype reference condition to be created
    PointVector mPointListDestination;                     /// A list that contents the all the points (from nodes) from the modelpart

    Properties::Pointer mpPairedProperties = nullptr;      /// This is the paired properties (unique for the given potential pair)

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method cleans the model part
     * @param rModelPart The model part of interest
     */
    virtual void CleanModelPart(ModelPart& rModelPart);

    /**
     * @brief This method checks the pairing
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     */
    virtual void CheckPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId
        );

    /**
     * @brief This method computes which nodes are active or inactive after after mapping the coordinates
     */
    virtual void ComputeActiveInactiveNodes();

    /**
     * @brief This method sets as active a node and it sets to an explicit approximation its LM
     * @param rNode The node reference to set
     * @param CommonEpsilon The penalty value
     * @param ScaleFactor The scale factor
     */
    virtual void SetActiveNode(
        Node& rNode,
        const double CommonEpsilon,
        const double ScaleFactor = 1.0
        );

    /**
     * @brief This method sets as inactive a node and it sets to zero its LM
     * @param ItNode The node reference to set
     */
    virtual void SetInactiveNode(Node& rNode);

    /**
     * @brief This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pObjectSlave The pointer to the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param pObjectMaster The pointer to the master condition
     * @param rMasterNormal The normal of the master condition
     * @param pIndexesPairs The map of indexes considered
     * @param pProperties The pointer to the Properties of the condition
     * @return The new created condition
     */
    virtual Condition::Pointer AddPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        GeometricalObject::Pointer pObjectSlave,
        const array_1d<double, 3>& rSlaveNormal,
        GeometricalObject::Pointer pObjectMaster,
        const array_1d<double, 3>& rMasterNormal,
        IndexMap::Pointer pIndexesPairs,
        Properties::Pointer pProperties
        );

    /**
     * @brief This converts the framework string to an enum
     * @param str The string
     * @return CheckGap: The equivalent enum
     */
    CheckGap ConvertCheckGap(const std::string& str);

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief This returns if we consider pure slip
     * @return True if we consider pure slip
     */
    bool IsPureSlip();

    /**
     * @brief This returns if we do not consider pure slip
     * @return True if we do not consider pure slip
     */
    bool IsNotPureSlip();

    /**
     * @brief This returns if we consider multiple searchs
     * @return True if we consider multiple searchs
     */
    bool IsMultipleSearchs();

    /**
     * @brief This returns if we do not consider multiple searchs
     * @return True if we do not consider multiple searchs
     */
    bool IsNotMultipleSearchs();

    /**
     * @brief This returns if we consider inverted search
     * @return True if we consider inverted search
     */
    bool IsInvertedSearch();

    /**
     * @brief This returns if we do not consider inverted search
     * @return True if we do not consider inverted search
     */
    bool IsNotInvertedSearch();

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

    /**
     * @brief This auxiliary method performs the seach using a KDTree
     * @param rSubContactModelPart The submodel part studied
     * @param rSubComputingContactModelPart The computing contact submodel part
     */
    void SearchUsingKDTree(
        ModelPart& rSubContactModelPart,
        ModelPart& rSubComputingContactModelPart
        );

    /**
     * @brief This auxiliary method performs the seach using a Octree
     * @param rSubContactModelPart The submodel part studied
     * @param rSubComputingContactModelPart The computing contact submodel part
     */
    void SearchUsingOcTree(
        ModelPart& rSubContactModelPart,
        ModelPart& rSubComputingContactModelPart
        );

    /**
     * @brief This method sets the origin destination model maps when only one model part is provided
     * @details The only model part should have MASTER/SLAVE flags in the nodes and conditions
     * @param rModelPart The main model part, where the origin/destination model parts will be created
     */
    void SetOriginDestinationModelParts(ModelPart& rModelPart);

    /**
     * @brief This function clears the mortar conditions already created
     * @param rNodesArray The array of nodes to clear
     */
    void ClearScalarMortarConditions(NodesArrayType& rNodesArray);

    /**
     * @brief This function clears the mortar conditions already created
     * @param rNodesArray The array of nodes to clear
     */
    void ClearComponentsMortarConditions(NodesArrayType& rNodesArray);

    /**
     * @brief This function clears the ALM frictionless mortar conditions already created
     * @param rNodesArray The array of nodes to clear
     */
    void ClearALMFrictionlessMortarConditions(NodesArrayType& rNodesArray);

    /**
     * @brief It check the conditions if they are correctly detected
     * @param pIndexesPairs Set containing the ids to the conditions
     * @param pGeometricalObject1 The pointer to the condition in the destination model part
     * @param pGeometricalObject2 The pointer to the condition in the destination model part
     * @param InvertedSearch If the search is inverted
     * @return If OK or Fail on the check
     */
    inline CheckResult CheckGeometricalObject(
        IndexMap::Pointer pIndexesPairs,
        const GeometricalObject::Pointer pGeometricalObject1,
        const GeometricalObject::Pointer pGeometricalObject2,
        const bool InvertedSearch = false
        );

    /**
     * @brief It check the conditions if they are correctly detected
     * @param pIndexesPairs Set containing the ids to the conditions
     * @param pCond1 The pointer to the condition in the destination model part
     * @param pCond2 The pointer to the condition in the destination model part
     * @param InvertedSearch If the search is inverted
     * @return If OK or Fail on the check
     */
    inline CheckResult CheckCondition(
        IndexMap::Pointer pIndexesPairs,
        const Condition::Pointer pCond1,
        const Condition::Pointer pCond2,
        const bool InvertedSearch = false
        );

    /**
     * @brief This method fills mPointListDestination
     */
    void FillPointListDestination();

    /**
     * @brief This method clears the destination list and
     * @param rSubContactModelPart The submodel part studied
     */
    void ClearDestinationListAndAssignFlags(ModelPart& rSubContactModelPart);

    /**
     * @brief This method computes search with KDTree
     * @param rTreePoints The tree points for search
     * @param rPointsFound The points found
     * @param rGeometry The geometry of the condition
     * @param TypeSearch The search type
     * @param SearchFactor The searh factor applied
     * @param AllocationSize The allocation size
     * @param Dynamic if the dynamic search is considered
     */
    inline IndexType PerformKDTreeSearch(
        KDTree& rTreePoints,
        PointVector& rPointsFound,
        GeometryType& rGeometry,
        const SearchTreeType TypeSearch = SearchTreeType::KdtreeInBox,
        const double SearchFactor = 3.5,
        const IndexType AllocationSize = 1000,
        const bool Dynamic = false
        );

    /**
     * @brief This method gets the maximum the ID of the conditions
     */
    inline IndexType GetMaximumConditionsIds();

    /**
     * @brief This method checks the potential pairing between two conditions/geometries (auxiliary one)
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pObjectSlave The pointer to the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param pObjectMaster The pointer to the master condition
     * @param rMasterNormal The normal of the master condition
     * @param pIndexesPairs The id sets of potential pairs
     * @param pProperties The pointer to the Properties of the condition
     * @param ActiveCheckFactor The value used auxiliarly to check if the node is in the potential contact zone
     * @param FrictionalProblem If the problem is frictional or not
     */
    void AddPotentialPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        GeometricalObject::Pointer pObjectSlave,
        const array_1d<double, 3>& rSlaveNormal,
        GeometricalObject::Pointer pObjectMaster,
        const array_1d<double, 3>& rMasterNormal,
        IndexMap::Pointer pIndexesPairs,
        Properties::Pointer pProperties,
        const double ActiveCheckFactor,
        const bool FrictionalProblem
        );

    /**
     * @brief This method computes the gap using a mapper
     * @param SearchOrientation The orientation of the search (inverted or not)
     */
    inline void ComputeMappedGap(const bool SearchOrientation);

    /**
     * @brief This method sets as inactive a node and it sets to zero its LM
     */
    inline void ComputeWeightedReaction();

    /**
     * @brief This method creates the auxiliary the pairing
     * @param rContactModelPart The modelpart  used in the assemble of the system
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     */
    inline void CreateAuxiliaryConditions(
        ModelPart& rContactModelPart,
        ModelPart& rComputingModelPart,
        IndexType& rConditionId
        );

    /**
     * @brief This method creates a debug file for normals
     * @param rModelPart The corresponding model part
     * @param rName The begining of the file name
     */
    void CreateDebugFile(
        ModelPart& rModelPart,
        const std::string& rName
        );

    /**
     * @brief Calculates the minimal distance between one node and its center
     * @return The radius of the geometry
     */
    static inline double Radius(GeometryType& ThisGeometry);

    /**
     * @brief This converts the framework string to an enum
     * @param str The string
     * @return SearchTreeType: The equivalent enum
     */
    SearchTreeType ConvertSearchTree(const std::string& str);

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

}; // Class BaseContactSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::istream& operator >> (std::istream& rIStream,
                                  BaseContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BaseContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.
