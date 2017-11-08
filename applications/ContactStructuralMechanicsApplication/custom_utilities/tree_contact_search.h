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
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

/* Custom includes*/
#include "custom_includes/point_item.h"

/* Custom utilities*/
#include "custom_utilities/search_utilities.h"

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
     * @param rMainModelPart: The model part to be considered
     * @param AllocationSize: The allocation considered in the search
     * @param ActiveCheckFactor: The factor considered to check if active or not
     * @param IntegrationOrder: The integration order considered
     * @param BucketSize: The size of the bucket
     * NOTE: Use an InterfacePreprocess object to create such a model part from a regular one:
     * InterfaceMapper = InterfacePreprocess()
     * InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    
    TreeContactSearch(
        ModelPart & rMainModelPart,
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
     * This function initializes the ALM frictionless mortar conditions already created 
     */
    
    void InitializeMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearScalarMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearComponentsMortarConditions();
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void TotalClearALMFrictionlessMortarConditions();
    
    /**
     * This function clears the mortar conditions already created (scalar version)
     */
    
    void PartialClearScalarMortarConditions();
    
    /**
     * This function clears the mortar conditions already created (components version)
     */
        
    void PartialClearComponentsMortarConditions();
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void PartialClearALMFrictionlessMortarConditions();
      
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
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param TypeSearch: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
    
    void UpdateMortarConditions();
    
    /**
     * This function has as pourpose to clean the existing pairs
     */
    
    void CleanMortarConditions();
    
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
    
    /**
     * This resets all the contact operators
     */
        
    void TotalResetContactOperators();
    
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
            
    /**
     * It check the conditions if they are correctly detected
     * @return ConditionPointers1: A vector containing the pointers to the conditions 
     * @param pCond1: The pointer to the condition in the destination model part
     * @param pCond2: The pointer to the condition in the destination model part  
     */
    
    static inline CheckResult CheckCondition(
        ConditionMap::Pointer& ConditionPointers1,
        const Condition::Pointer& pCond1,
        const Condition::Pointer& pCond2,
        const bool InvertedSearch = false
        );
    
    /**  
     * Calculates the minimal distance between one node and its center 
     * @return The radius of the geometry 
     */ 
    
    static inline double Radius(GeometryType& ThisGeometry);
    
    /**
     * This converts the framework string to an enum
     * @param str: The string
     * @return SearchTreeType: The equivalent enum
     */
    
    SearchTreeType ConvertSearchTree(const std::string& str);
    
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
  
    ModelPart& mrMainModelPart;               // The main model part
    unsigned int mDimension;                  // Dimension size of the space
    unsigned int mAllocationSize;             // Allocation size for the vectors and max number of potential results
    double mSearchFactor;                     // The search factor to be considered
    bool mDualSearchCheck;                    // The search is done reciprocally
    bool mStrictSearchCheck;                  // The search is done requiring IsInside as true
    bool mUseExactIntegration;                // The search filter the results with the exact integration
    bool mInvertedSearch;                     // The search will be done inverting the way master and slave/master is assigned
    SearchTreeType mSearchTreeType;           // The search tree considered
    unsigned int mBucketSize;                 // Bucket size for kd-tree
    PointVector mPointListDestination;        // A list that contents the all the points (from nodes) from the modelpart 

    ///@}
    ///@name Private Operators
    ///@{

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

    ///@}

}; // Class TreeContactSearch

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  TreeContactSearch& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TreeContactSearch& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
