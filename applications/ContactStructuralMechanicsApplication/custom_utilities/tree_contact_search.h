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
    
///@}
///@name  Enum's
///@{
    
    enum SearchTreeType {KdtreeInRadius = 0, KdtreeInBox = 1, Kdop = 2};
    
    enum CheckResult {Fail = 0, AlreadyInTheMap = 1, OK = 2};
    
    enum CheckGap {NoCheck = 0, DirectCheck = 1, MappingCheck = 2};
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** \brief TreeContactSearch
 * This utilitiy has as objective to create the contact conditions.
 * The conditions that can be created are Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 */
template<unsigned int TDim, unsigned int TNumNodes>
class TreeContactSearch
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    
    // Type definitions for the tree
    typedef PointItem                                             PointType;
    typedef PointType::Pointer                             PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of TreeContactSearch
    KRATOS_CLASS_POINTER_DEFINITION( TreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    
    /**
     * The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part to be considered
     * @param ThisParameters The condiguration parameters, it includes:
     *                       - The allocation considered in the search
     *                       - The factor considered to check if active or not
     *                       - The integration order considered
     *                       - The size of the bucket
     *                       - The proportion increased of the Radius/Bounding-box volume for the search
     *                       - TypeSearch: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * NOTE: Use an InterfacePreprocess object to create such a model part from a regular one:
     * InterfaceMapper = InterfacePreprocess()
     * InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    
    TreeContactSearch( ModelPart& rMainModelPart, Parameters ThisParameters =  Parameters(R"({})") ); 
    
    virtual ~TreeContactSearch()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function initializes the ALM frictionless mortar conditions already created 
     */
    
    void InitializeMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void ClearScalarMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void ClearComponentsMortarConditions();
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void ClearALMFrictionlessMortarConditions();
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar();

    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar();

    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     */
    
    void UpdateMortarConditions();
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions();
    
    /**
     * It sets if the search is inverted
     */
    
    void InvertSearch();
    
    /**
     * This resets the contact operators
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
  
    ModelPart& mrMainModelPart;                      // The main model part
    Parameters mThisParameters;                      // The configuration parameters
    CheckGap mCheckGap;                              // If the gap is checked during the search
    bool mInvertedSearch;                            // The search will be done inverting the way master and slave/master is assigned
    std::string mConditionName;                      // The name of the condition to be created
    bool mCreateAuxiliarConditions;                  // If the auxiliar conditions are created or not
    PointVector mPointListDestination;               // A list that contents the all the points (from nodes) from the modelpart 

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
       
    /**
     * This method computes the maximal nodal H
     */
    inline double GetMaxNodalH();
       
    /**
     * This method computes the mean nodal H
     */
    inline double GetMeanNodalH();
       
    /**
     * This method initializes the acceleraction when there is a volume acceleration
     */
    inline void InitializeAcceleration();
    
    /**
     * It check the conditions if they are correctly detected
     * @return ConditionPointers1: A vector containing the pointers to the conditions 
     * @param pCond1 The pointer to the condition in the destination model part
     * @param pCond2 The pointer to the condition in the destination model part  
     * @param InvertedSearch If the search is inverted
     */
    
    static inline CheckResult CheckCondition(
        IndexSet::Pointer IndexesSet,
        const Condition::Pointer pCond1,
        const Condition::Pointer pCond2,
        const bool InvertedSearch = false
        );
    
    /**
     * This method reorders the ID of the conditions
     */

    inline std::size_t ReorderConditionsIds();
    
    /**
     * This method checks the potential pairing between two conditions/geometries
     */
    inline void AddPotentialPairing(
        ModelPart& rComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        PointVector& PointsFound,
        const unsigned int NumberOfPointsFound,
        IndexSet::Pointer IndexesSet
        );
    
    /**
     * This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param ConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param pCondMaster The pointer to the master condition
     */
    inline void AddPairing(
        ModelPart& rComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        Condition::Pointer pCondMaster
        );
    
    /**
     * This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param ConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param pCondMaster The pointer to the master condition
     * @param IndexesSet The map of indexes considered
     */
    inline void AddPairing(
        ModelPart& rComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        Condition::Pointer pCondMaster,
        IndexSet::Pointer IndexesSet
        );
    
    /**
     * This method checks the pairing
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param ConditionId The ID of the new condition to be created
     */
    inline void CheckPairing(
        ModelPart& rComputingModelPart,
        std::size_t& ConditionId
        );
    
    /**
     * This method creates the auxiliar the pairing
     * @param rContactModelPart The modelpart  used in the assemble of the system
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param ConditionId The ID of the new condition to be created
     */
    inline void CreateAuxiliarConditions(
        ModelPart& rContactModelPart,
        ModelPart& rComputingModelPart,
        std::size_t& ConditionId
        );
    
    /**  
     * Calculates the minimal distance between one node and its center 
     * @return The radius of the geometry 
     */ 
    
    static inline double Radius(GeometryType& ThisGeometry);
    
    /**
     * This converts the framework string to an enum
     * @param str The string
     * @return SearchTreeType: The equivalent enum
     */
    
    SearchTreeType ConvertSearchTree(const std::string& str);
    
    /**
     * This converts the framework string to an enum
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
