// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/search_utilities.h"
#include "custom_utilities/point_item.h"

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
        )
    :mrMainModelPart(rMainModelPart.GetSubModelPart("Contact")),
     mDimension(rMainModelPart.GetProcessInfo()[DOMAIN_SIZE])
    {        
        Parameters DefaultParameters = Parameters(R"(
        {
            "allocation_size"                      : 1000, 
            "bucket_size"                          : 4, 
            "search_factor"                        : 2.0, 
            "type_search"                          : "InRadius", 
            "active_check_factor"                  : 0.01,
            "dual_search_check"                    : false,
            "strict_search_check"                  : true
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mAllocationSize = ThisParameters["allocation_size"].GetInt();
        mSearchFactor = ThisParameters["search_factor"].GetDouble();
        mActiveCheckFactor = ThisParameters["active_check_factor"].GetDouble();
        mDualSearchCheck = ThisParameters["dual_search_check"].GetBool();
        mStrictSearchCheck = ThisParameters["strict_search_check"].GetBool();
        mSearchTreeType = ConvertSearchTree(ThisParameters["type_search"].GetString());
        mBucketSize = ThisParameters["bucket_size"].GetInt();
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->Set(ACTIVE, false);
        }
        
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

        #pragma omp parallel for 
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            itCond->Set(ACTIVE, false);
        }
    }
    
    virtual ~TreeContactSearch(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function initializes the ALM frictionless mortar conditions already created 
     */
    
    void InitializeMortarConditions()
    {
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

//         #pragma omp parallel for 
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;

            itCond->GetValue(CONTACT_SETS) = boost::shared_ptr<ConditionSet>(new ConditionSet); 
//             itCond->GetValue(CONTACT_SETS)->reserve(mAllocationSize); 
        }
    }
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearScalarMortarConditions()
    {
        ResetContactOperators(mrMainModelPart);
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {
                itNode->Set( ACTIVE, false );
                itNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
            }
        }  
    }
    
    void TotalClearComponentsMortarConditions()
    {
        ResetContactOperators(mrMainModelPart);
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {
                itNode->Set( ACTIVE, false );
                itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER) = ZeroVector(3);
            }
        }  
    }
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void TotalClearALMFrictionlessMortarConditions()
    {
        ResetContactOperators(mrMainModelPart);
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {
                itNode->Set( ACTIVE, false );
                itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
            }
        }  
    }
    
    /**
     * This function clears the mortar conditions already created (scalar version)
     */
    
    void PartialClearScalarMortarConditions()
    {
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
            }
        } 
    }
    
    /**
     * This function clears the mortar conditions already created (components version)
     */
        
    void PartialClearComponentsMortarConditions()
    {
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER) = ZeroVector(3);
            }
        } 
    }
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void PartialClearALMFrictionlessMortarConditions()
    {
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
            }
        } 
    }
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar()
    {
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

        #pragma omp for nowait schedule(static)
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            if (itCond->Is(MASTER) == true)
            {
                PointTypePointer pPoint = PointTypePointer(new PointItem((*itCond.base())));
                (mPointListDestination).push_back(pPoint);
            }
        }
    }

    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar()
    {
        auto numPoints = mPointListDestination.end() - mPointListDestination.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numPoints; i++) 
        {
            mPointListDestination[i]->UpdatePoint();
        }

    }

    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param TypeSearch: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
    
    void UpdateMortarConditions()
    {        
        // We update the list of points
        UpdatePointListMortar();
        
        // Calculate the mean of the normal in all the nodes
        ContactUtilities::ComputeNodesMeanNormalModelPart(mrMainModelPart); 
        
//         #pragma omp parallel 
//         {
            // Initialize values
            PointVector PointsFound(mAllocationSize);
            unsigned int NumberPointsFound = 0;    
            
            // Create a tree
            // It will use a copy of mNodeList (a std::vector which contains pointers)
            // Copying the list is required because the tree will reorder it for efficiency
            KDTree TreePoints(mPointListDestination.begin(), mPointListDestination.end(), mBucketSize);
            
            // Iterate in the conditions
            ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
            auto numConditions = pConditions.end() - pConditions.begin();

//             #pragma omp for 
            for(unsigned int i = 0; i < numConditions; i++) 
            {
                auto itCond = pConditions.begin() + i;
                
                if (itCond->Is(SLAVE) == true)
                {
                    if (mSearchTreeType == KdtreeInRadius)
                    {
                        Point<3> Center;
                        const double SearchRadius = mSearchFactor * ContactUtilities::CenterAndRadius((*itCond.base()), Center);

                        NumberPointsFound = TreePoints.SearchInRadius(Center, SearchRadius, PointsFound.begin(), mAllocationSize);
                    }
                    else if (mSearchTreeType == KdtreeInBox)
                    {
                        const Point<3> Center = itCond->GetGeometry().Center();
                        Node<3> MinPoint, MaxPoint;
                        itCond->GetGeometry().BoundingBox(MinPoint, MaxPoint);
                        ContactUtilities::ScaleNode<Node<3>>(MinPoint, Center, mSearchFactor);
                        ContactUtilities::ScaleNode<Node<3>>(MaxPoint, Center, mSearchFactor);
                        NumberPointsFound = TreePoints.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), mAllocationSize);
                    }
                    else
                    {
                        KRATOS_ERROR << " The type search is not implemented yet does not exist!!!!. SearchTreeType = " << mSearchTreeType << std::endl;
                    }
                    
                    if (NumberPointsFound > 0)
                    {                           
                        boost::shared_ptr<ConditionSet>& ConditionPointersDestination = itCond->GetValue(CONTACT_SETS);
                        
                        for(unsigned int i = 0; i < NumberPointsFound; i++)
                        {   
                            Condition::Pointer pCondOrigin = PointsFound[i]->GetCondition();
                            
                            if (CheckCondition(ConditionPointersDestination, (*itCond.base()), pCondOrigin) == true) 
                            {    
                                // If not active we check if can be potentially in contact
                                SearchUtilities::ContactContainerFiller(ConditionPointersDestination, (*itCond.base()), pCondOrigin, itCond->GetValue(NORMAL), pCondOrigin->GetValue(NORMAL), mActiveCheckFactor, mDualSearchCheck, mStrictSearchCheck); 
                            }
                            
                            if (ConditionPointersDestination->size() > 0)
                            {                        
                                itCond->Set(ACTIVE, true);
                            }
                        }
                    }
                }
            }
//         }
    }
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions()
    {
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

//         #pragma omp parallel for 
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            if (itCond->Is(MASTER) == true && itCond->Is(ACTIVE) == true)
            {
                KRATOS_WATCH(itCond->Id());
                KRATOS_ERROR << "THIS IS NOT SUPPOSED TO HAPPEN" << std::endl; 
            }
            
            if (itCond->Is(SLAVE) == true && itCond->Is(ACTIVE) == true)
            {
                KRATOS_WATCH(itCond->GetGeometry());
                
                boost::shared_ptr<ConditionSet>& ConditionPointersDestination = itCond->GetValue(CONTACT_SETS);
                KRATOS_WATCH(ConditionPointersDestination->size());
                ConditionPointersDestination->print();
            }
        }
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
//         #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {
                std::cout << "Node: " << itNode->Id() << " is active" << std::endl;
            }
        }
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
    
    bool CheckCondition(
        boost::shared_ptr<ConditionSet>& ConditionPointers1,
        const Condition::Pointer & pCond1,
        const Condition::Pointer & pCond2
        )
    {
        if (((pCond1 != pCond2) && (pCond1->GetValue(ELEMENT_POINTER) != pCond2->GetValue(ELEMENT_POINTER))) == false) // Avoiding "auto self-contact" and "auto element contact"
        {
            return false;
        }
        
        // Avoid conditions oriented in the same direction
        const double Tolerance = 1.0e-16;
        if (norm_2(pCond1->GetValue(NORMAL) - pCond2->GetValue(NORMAL)) < Tolerance)
        {
            return false;
        }
        
        if (ConditionPointers1->find(pCond2) != ConditionPointers1->end())
        {
            return false;
        }
        
        if (pCond2->Is(SLAVE) == true) // Otherwise will not be necessary to check
        {
            auto& ConditionPointers2 = pCond2->GetValue(CONTACT_SETS);
            
            if (ConditionPointers2->find(pCond1) != ConditionPointers2->end())
            {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * This resets the contact operators
     * @param rModelPart: The model part where the contact operators are reset
     */
        
    void ResetContactOperators(ModelPart & rModelPart)
    {
        ConditionsArrayType& pCond = rModelPart.Conditions();
    
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            if (itCond->Is(SLAVE) == true && itCond->Is(ACTIVE) == true)
            {
                itCond->Set(ACTIVE, false);
                
                auto& ConditionPointers = itCond->GetValue(CONTACT_SETS);
                
                if (ConditionPointers != nullptr)
                {
                    ConditionPointers->clear();
//                     ConditionPointers->reserve(mAllocationSize); 
                }
            }
        }   
    }
    
    /**
     * This converts the framework string to an enum
     * @param str: The string
     * @return SearchTreeType: The equivalent enum
     */
    
    SearchTreeType ConvertSearchTree(const std::string& str)
    {
        if(str == "InRadius") 
        {
            return KdtreeInRadius;
        }
        else if(str == "InBox") 
        {
            return KdtreeInBox;
        }
        else if (str == "KDOP")
        {
            KRATOS_ERROR << "KDOP contact search: Not yet implemented" << std::endl;
            return Kdop;
        }
        else
        {
            return KdtreeInRadius;
        }
    }
    
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
    double mActiveCheckFactor;                // The check factor to be considered
    bool mDualSearchCheck;                    // The search is done reciprocally
    bool mStrictSearchCheck;                  // The search is done requiring IsInside as true
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
