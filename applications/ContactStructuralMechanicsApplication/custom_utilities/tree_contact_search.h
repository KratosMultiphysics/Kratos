// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/simple_mortar_mapper_process.h" 
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

/* Custom includes*/
#include "custom_includes/point_item.h"
#include "custom_conditions/paired_condition.h"

/* Tree structures */
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    typedef std::size_t SizeType;

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
 * @class TreeContactSearch 
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This utilitiy has as objective to create the contact conditions.
 * @details The conditions that can be created are Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
class TreeContactSearch
{
public:
    ///@name Type Definitions
    ///@{
    
    /// General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    
    /// Index type definition
    typedef std::size_t                                           IndexType;

    /// Type definitions for the tree
    typedef PointItem                                             PointType;
    typedef PointType::Pointer                             PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;
    
    /// KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// The type of mapper considered
    typedef SimpleMortarMapperProcess<TDim, TNumNodes, Variable<array_1d<double, 3>>, TNumNodesMaster> MapperType;

    /// The definition of zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of TreeContactSearch
    KRATOS_CLASS_POINTER_DEFINITION( TreeContactSearch );
      
    ///@}
    ///@name  Enum's
    ///@{

    enum class SearchTreeType {KdtreeInRadius = 0, KdtreeInBox = 1, Kdop = 2};

    enum class CheckResult {Fail = 0, AlreadyInTheMap = 1, OK = 2};

    enum class CheckGap {NoCheck = 0, DirectCheck = 1, MappingCheck = 2};

    enum class TypeSolution {NormalContactStress = 0, ScalarLagrangeMultiplier = 1, VectorLagrangeMultiplier = 2};

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
     * @todo Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @note Use an InterfacePreprocess object to create such a model part from a regular one:
     *          -# InterfaceMapper = InterfacePreprocess()
     *          -# InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    TreeContactSearch( 
        ModelPart& rMainModelPart, 
        Parameters ThisParameters =  Parameters(R"({})") 
        ); 
    
    virtual ~TreeContactSearch()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function initializes the ALM frictionless mortar conditions already created 
     */
    void InitializeMortarConditions();
    
    /**
     * @brief This function clears the mortar conditions already created 
     */
    void ClearMortarConditions();
      
    /**
     * @brief This method checks that the contact model part is unique (so the model parts contain unique contact pairs)
     */
    void CheckContactModelParts();

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
    void ResetContactOperators();
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    virtual std::string Info() const
    {
        return "TreeContactSearch";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    virtual void PrintInfo(std::ostream& rOStream) const
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
  
    ModelPart& mrMainModelPart;        /// The main model part
    Parameters mThisParameters;        /// The configuration parameters
    CheckGap mCheckGap;                /// If the gap is checked during the search
    TypeSolution mTypeSolution;        /// The solution type 
    bool mInvertedSearch;              /// The search will be done inverting the way master and slave/master is assigned
    std::string mConditionName;        /// The name of the condition to be created
    bool mCreateAuxiliarConditions;    /// If the auxiliar conditions are created or not
    PointVector mPointListDestination; /// A list that contents the all the points (from nodes) from the modelpart 
    bool mMultipleSearchs;             /// If we consider multiple serach or not
    bool mPredefinedMasterSlave;       /// If the master/slave sides are predefined

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
       
    /**
     * @brief This method sets the origin destination model maps when only one model part is provided
     * @details The only model part should have MASTER/SLAVE flags in the nodes and conditions
     * @param rModelPart The main model part, where the origin/destination model parts will be created
     */
    void SetOriginDestinationModelParts(ModelPart& rModelPart);

    /**
     * @brief This function clears the mortar conditions already created 
     * @param NodesArray The array of nodes to clear
     */
    void ClearScalarMortarConditions(NodesArrayType& NodesArray);
    
    /**
     * @brief This function clears the mortar conditions already created 
     * @param NodesArray The array of nodes to clear
     */
    void ClearComponentsMortarConditions(NodesArrayType& NodesArray);
    
    /**
     * @brief This function clears the ALM frictionless mortar conditions already created 
     * @param NodesArray The array of nodes to clear
     */
    void ClearALMFrictionlessMortarConditions(NodesArrayType& NodesArray);
       
    /**
     * @brief This computes a simple linear regression to the gap and contact pressure
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void ComputeLinearRegressionGapPressure(
        double& a,
        double& b
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
     * @brief This method is used in case of not predefined master/slave we assign the master/slave nodes and conditions
     * @param rModelPart The model part to assign the flags
     */
    static inline void NotPredefinedMasterSlave(ModelPart& rModelPart);

    /**
     * @brief This method gets the maximum the ID of the conditions
     */
    inline IndexType GetMaximumConditionsIds();
    
    /**
     * @brief This method checks the potential pairing between two conditions/geometries
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param rPointsFound The potential pairs found 
     * @param NumberOfPointsFound The number of potential pairs found
     * @param IndexesPairs The id sets of potential pairs
     */
    inline void AddPotentialPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        Condition::Pointer pCondSlave,
        PointVector& rPointsFound,
        const IndexType NumberOfPointsFound,
        IndexMap::Pointer IndexesPairs
        );
    
    /**
     * @brief This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param pCondMaster The pointer to the master condition
     */
    inline void AddPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        Condition::Pointer pCondSlave,
        Condition::Pointer pCondMaster
        );
    
    /**
     * @brief This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param pCondMaster The pointer to the master condition
     * @param IndexesPairs The map of indexes considered
     */
    inline void AddPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        Condition::Pointer pCondSlave,
        Condition::Pointer pCondMaster,
        IndexMap::Pointer IndexesPairs
        );
    
    /**
     * @brief This method checks the pairing
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     */
    inline void CheckPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId
        );
    
    /**
     * @brief This method computes the gap using a mapper 
     * @param SearchOrientation The orientation of the search (inverted or not)
     */
    inline void ComputeMappedGap(const bool SearchOrientation);
    
    /**
     * @brief This method computes which nodes are active or inactive after after mapping the coordinates
     */
    inline void ComputeActiveInactiveNodes();
    
    /**
     * @brief This method sets as active a node and it predicts the value of its LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void SetActiveNode(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This method sets as inactive a node and it sets to zero its LM
     * @param ItNode The node iterator to set
     */
    inline void SetInactiveNode(NodesArrayType::iterator ItNode);
    
    /**
     * @brief This function predicts the scalar LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void CorrectScalarMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the vector LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void CorrectComponentsMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictionless LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void CorrectALMFrictionlessMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictionless in components LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void CorrectALMFrictionlessComponentsMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictional LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void CorrectALMFrictionalMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the scalar LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void PredictScalarMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the vector LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void PredictComponentsMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictionless LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void PredictALMFrictionlessMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictionless in components LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void PredictALMFrictionlessComponentsMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This function predicts the ALM frictional LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    inline void PredictALMFrictionalMortarLM(
        NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );
    
    /**
     * @brief This method sets as inactive a node and it sets to zero its LM
     */
    inline void ComputeWeightedReaction();
    
    /**
     * @brief This method switchs the flag of an array of nodes
     * @param rNodes The set of nodes where the flags are reset
     */
    static inline void SwitchFlagNodes(NodesArrayType& rNodes)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rNodes.size()); ++i) {
            auto it_node = rNodes.begin() + i;
            it_node->Flip(SLAVE);
            it_node->Flip(MASTER);
        }
    }
    
    /**
     * @brief This method creates the auxiliar the pairing
     * @param rContactModelPart The modelpart  used in the assemble of the system
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     */
    inline void CreateAuxiliarConditions(
        ModelPart& rContactModelPart,
        ModelPart& rComputingModelPart,
        IndexType& rConditionId
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
    
    /**
     * @brief This converts the framework string to an enum
     * @param str The string
     * @return CheckGap: The equivalent enum
     */
    CheckGap ConvertCheckGap(const std::string& str);
    
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

}; // Class TreeContactSearch

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

// /****************************** INPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   TreeContactSearch& rThis);
// 
// /***************************** OUTPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const TreeContactSearch& rThis)
// {
//     return rOStream;
// }

///@}

}  // namespace Kratos.

#endif // KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
