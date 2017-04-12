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

#if !defined(KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "includes/model_part.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/search_utilities.h"
#include "custom_utilities/point_item.h"

// TODO: Add parallelization

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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief DeprecatedTreeContactSearch NOTE: OLD VERSION
 * This utilitiy has as objective to create the contact conditions. The conditions that can be created are:
 * * Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments (En principio)
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree (En principio no)
 * To consider autocontact use the same model_part as origin and destination (En principio, preguntar a Riccardo)
 */

class DeprecatedTreeContactSearch
{

public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType               NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef GeometryData::IntegrationMethod             IntegrationMethod;
    typedef Node<3>                                     NodeType;
    typedef Geometry<NodeType>                          GeometryType;
    
    // Type definitions for the tree
    typedef PointItem                                    PointType;
    typedef PointItem::Pointer                           PointTypePointer;
    typedef std::vector<PointType::Pointer>              PointVector;
    typedef std::vector<PointType::Pointer>::iterator    PointIterator;
    typedef std::vector<double>                          DistanceVector;
    typedef std::vector<double>::iterator                DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > tree;

    /// Pointer definition of DeprecatedTreeContactSearch
    // KRATOS_CLASS_POINTER_DEFINITION( DeprecatedTreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{
    
    // Class Constructor
    // WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
    // Use an InterfacePreprocess object to create such a model part from a regular one:
    // InterfaceMapper = InterfacePreprocess()
    // InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
    DeprecatedTreeContactSearch(
        ModelPart & rOriginModelPart,
        ModelPart & rDestinationModelPart,
        const unsigned int allocation_size
        ):
        mrOriginModelPart(rOriginModelPart),
        mrDestinationModelPart(rDestinationModelPart),
        mBucketSize(4),
        mdimension(rOriginModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension()),
        mallocation(allocation_size)
    {  
    //     // Destination model part
    //     ModelPartSetter(mrDestinationModelPart, false, true, false);
        
        NodesArrayType& pNode = mrDestinationModelPart.Nodes();
            
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->Set(SLAVE, true);  
        }

    //     // Origin model part
    //     ModelPartSetter(mrOriginModelPart, false, false, true); 
        
        ConditionsArrayType& pCond = mrOriginModelPart.Conditions();
            
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            itCond->Set(MASTER, true);  
        }
        
    }
    
    // Class Destructor
    virtual ~DeprecatedTreeContactSearch(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function sets the flags in the model part
     */
    void ModelPartSetter(
        ModelPart& rModelPart,
        const bool rActive,
        const bool rSlave,
        const bool rMaster
        ) 
    {
        NodesArrayType& pNode = rModelPart.Nodes();
            
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->Set( SLAVE,   rSlave );  
            itNode->Set( ACTIVE, rActive );  // NOTE: It is supposed to be already false, just in case  
        }
        
        ConditionsArrayType& pCond  = rModelPart.Conditions();
            
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            itCond->Set( ACTIVE, rActive ); // NOTE: It is supposed to be already false, just in case   
            itCond->Set( MASTER, rMaster);
        }
    }
    
    /**
     * This function initializes the mortar conditions already created 
     */
    
    void InitializeMortarConditions(
        const double rActiveCheckFactor,
        const double rAugmentationNormal,
        const double rAugmentationTangent,
        const int rIntegrationOrder
        )
    {
        // Destination model part
        InitializeConditions(mrDestinationModelPart, rActiveCheckFactor, rAugmentationNormal, rAugmentationTangent, rIntegrationOrder);
        
        // Origin model part
        InitializeConditions(mrOriginModelPart, rActiveCheckFactor, rAugmentationNormal, rAugmentationTangent, rIntegrationOrder);
    }

    /**
     * This function initializes conditions
     */
    
    void InitializeConditions(
        ModelPart & rModelPart, 
        const double rActiveCheckFactor,
        const double rAugmentationNormal,
        const double rAugmentationTangent,
        const int rIntegrationOrder
        )
    {   
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCond.ptr_end();

        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
//             cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
            cond_it->GetProperties().SetValue(ACTIVE_CHECK_FACTOR, rActiveCheckFactor);
            cond_it->GetProperties().SetValue(NORMAL_AUGMENTATION_FACTOR,  rAugmentationNormal);
            cond_it->GetProperties().SetValue(TANGENT_AUGMENTATION_FACTOR, rAugmentationTangent);
            if (cond_it->GetProperties().Has(INTEGRATION_ORDER_CONTACT) == false)
            {
                cond_it->GetProperties().SetValue(INTEGRATION_ORDER_CONTACT, rIntegrationOrder);
            }
        }
    }
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearMortarConditions()
    {
        // Destination model part
        TotalClearConditions(mrDestinationModelPart);
    }

    /**
     * This function clears the mortar conditions already created 
     */
    
    void PartialClearMortarConditions()
    {
        // Destination model part
        PartialClearConditions(mrDestinationModelPart);
    }
    
    /**
     * This function clears conditions already created 
     */
    
    void TotalClearConditions(ModelPart & rModelPart)
    {
        ResetContactOperators(rModelPart);
        
        NodesArrayType& pNode = rModelPart.Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == true)
            {
                itNode->Set( ACTIVE, false );
                itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
            }
        }
    }
    
    /**
     * This function clears partially the conditions already created 
     */
    
    void PartialClearConditions(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
            }
        }
    }
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar()
    {
        // Destination model part
        CreatePointListConditions(mrDestinationModelPart, mPointListDestination);
    }

    /**
     * This function creates a condition list
     */
    
    void CreatePointListConditions(
        ModelPart & rModelPart, 
        PointVector & PoinList
        )
    {
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCond.ptr_end();
        
        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            const Condition::Pointer & pCond = (*cond_it.base());
            ContactUtilities::ConditionNormal(pCond);
            PointItem::Pointer pPoint = PointItem::Pointer(new PointItem(pCond));
            (PoinList).push_back(pPoint);
        }
    }
    
    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar()
    {
        // Destination model part
        UpdatePointListConditions(mrDestinationModelPart, mPointListDestination);
    }

    /**
     * This function updates a condition list
     */
    
    void UpdatePointListConditions(
        ModelPart & rModelPart, 
        PointVector & PoinList
        )
    {
        auto numPoints = PoinList.end() - PoinList.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numPoints; i++) 
        {
            auto itPoint = PoinList.begin() + i;
            
            (*itPoint.base())->UpdatePoint();
        }
    }
    
    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
        
    void CreateMortarConditions(
        const double SearchFactor,
        const int type_search
    )
    {
        TotalClearMortarConditions(); // Clear the conditions
        UpdateMortarConditions(SearchFactor, type_search);
    }
    
    void UpdateMortarConditions(
        const double SearchFactor,
        const int type_search
    )
    {
        // Initialize values
        PointVector PointsFound(mallocation);
        std::vector<double> PointsDistances(mallocation);
        unsigned int NumberPointsFound = 0;    
        
        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        tree Tree_points(mPointListDestination.begin(), mPointListDestination.end(), mBucketSize);
        
        ConditionsArrayType& pCond  = mrOriginModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCond.ptr_end();
        
        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            const Condition::Pointer pCondOrigin = (*cond_it.base());
            
            if (type_search == 0)
            {
                Point<3> Center;
                const double SearchRadius = SearchFactor * ContactUtilities::CenterAndRadius(pCondOrigin, Center);
                ContactUtilities::ConditionNormal(pCondOrigin);

                NumberPointsFound = Tree_points.SearchInRadius(Center, SearchRadius, PointsFound.begin(), PointsDistances.begin(), mallocation);
            }
            else if (type_search == 1)
            {
                Point<3> Center;
                ContactUtilities::CenterAndRadius(pCondOrigin, Center);
                ContactUtilities::ConditionNormal(pCondOrigin);
                
                Node<3> MinPoint, MaxPoint;
                ContactUtilities::ScaleNode<Node<3>>(MinPoint, Center, SearchFactor);
                ContactUtilities::ScaleNode<Node<3>>(MaxPoint, Center, SearchFactor);
                pCondOrigin->GetGeometry().BoundingBox(MinPoint, MaxPoint);
                NumberPointsFound = Tree_points.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), mallocation);
            }
            else
            {
                KRATOS_ERROR << " The type search declared does not exist!!!!. type_search = " << type_search << std::endl;
            }
            
            if (NumberPointsFound > 0)
            {
//                 KRATOS_WATCH(NumberPointsFound); 
                for(unsigned int i = 0; i < NumberPointsFound; i++)
                {   
                    Condition::Pointer pCondDestination = PointsFound[i]->GetCondition();
                    
                    std::vector<contact_container> *& ConditionPointersDestination = pCondDestination->GetValue(CONTACT_CONTAINERS);
                    
                    bool to_check_cond = CheckCondition(ConditionPointersDestination, pCondDestination, pCondOrigin);
                    
                    if (to_check_cond == true) 
                    {    
                        // If not active we check if can be potentially in contact
                        MortarContainerFiller(pCondDestination, pCondOrigin, ConditionPointersDestination, pCondDestination->GetProperties().GetValue(ACTIVE_CHECK_FACTOR));
                    }
                }
            }
        }
        
        // Here we remove all the inactive pairs
        ClearAllInactivePairs(mrDestinationModelPart); 
        
        // Calculate the mean of the normal in all the nodes
        ComputeNodesMeanNormal();
    }
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions()
    {    
        ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCondDestination.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCondDestination.ptr_end();
        
        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            if (cond_it->Is(ACTIVE))
            {
                std::vector<contact_container> *& ConditionPointersDestination = cond_it->GetValue(CONTACT_CONTAINERS);
                KRATOS_WATCH(ConditionPointersDestination->size());
                
                for (unsigned int i = 0; i < ConditionPointersDestination->size(); i++)
                {
                    (*ConditionPointersDestination)[i].print();
                } 
            }
        }
        
        NodesArrayType& pNodeDestination    = mrDestinationModelPart.Nodes();
        NodesArrayType::iterator node_begin = pNodeDestination.ptr_begin();
        NodesArrayType::iterator node_end   = pNodeDestination.ptr_end();
        
        for(NodesArrayType::iterator node_it = node_begin; node_it!=node_end; node_it++)
        {         
            if (node_it->Is(ACTIVE) == true)
            {
                std::cout << "Node: " << node_it->Id() << " is active" << std::endl;
            }
        }
    }
    
    /**
     * It clears all the inactive pairs
     * @param rModelPart: The modelpart to clear
     */
    
    void ClearAllInactivePairs(ModelPart & rModelPart)
    {
        ConditionsArrayType& pCond = rModelPart.Conditions();
        
        auto numConditions = pCond.end() - pCond.begin();
        
    //     #pragma omp parallel for // NOTE: Be careful, if you change something with get value over the nodes iteraring in the conditions this will not work in OpenMP 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            if (itCond->Is(ACTIVE) == true)
            {            
                std::vector<contact_container> *& ConditionPointers = itCond->GetValue(CONTACT_CONTAINERS);
                                
                if (ConditionPointers != NULL)
                {
                    std::vector<contact_container> AuxConditionPointers;
                    for (unsigned int pair = 0; pair < (*ConditionPointers).size();pair++)
                    {
                        contact_container aux_contact_container = (*ConditionPointers)[pair];
                        
                        if (aux_contact_container.active_pair == false)
                        {
                            // Last oportunity for the condition pair
                            const double ActiveCheckFactor = itCond->GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
                            const bool cond_active = SearchUtilities::ContactChecker(itCond->GetGeometry(),    (aux_contact_container.condition)->GetGeometry(), 
                                                                                            itCond->GetValue(NORMAL), (aux_contact_container.condition)->GetValue(NORMAL), 
                                                                                            ActiveCheckFactor);
                            if (cond_active == true) // Still paired
                            {
                                aux_contact_container.active_pair = true; 
                                AuxConditionPointers.push_back(aux_contact_container);
                            }
                        }
                        else // It is already active pair, we append
                        {
                            AuxConditionPointers.push_back(aux_contact_container);
                        }
                    } 
                    
                    // Now we copy the final result
                    *ConditionPointers = AuxConditionPointers;
                    
                    // All the pairs has been removed
                    if ((*ConditionPointers).size() == 0)
                    {
                        itCond->Set(ACTIVE, false);
                    }
                }
                else
                {
                    itCond->Set(ACTIVE, false);
                }

            }
        }
    }
    
    /**
     * It check the conditions if they are correctly detected
     * @return ConditionPointers: A vector containing the pointers to the conditions 
     * @param pCondDestination: The pointer to the condition in the destination model part
     * @param pCondOrigin: The pointer to the condition in the destination model part  
     */
    
    bool CheckCondition(
        std::vector<contact_container> *& ConditionPointers,
        const Condition::Pointer & pCondDestination,
        const Condition::Pointer & pCondOrigin
        )
    {
        bool aux_bool = (pCondDestination != pCondOrigin); // Avoiding "auto self-contact"
        
        if (aux_bool == true)
        {
            for (unsigned int pair = 0; pair < ConditionPointers->size(); pair++)
            {
                if ((*ConditionPointers)[pair].condition == pCondOrigin)
                {
                    aux_bool = false;
                    break;
                }
            }
        }
        
        return aux_bool;
    }
    
    /**
     * Fills the contact container variable
     * @param pCond_1: The origin condition
     * @param pCond_2: The potential condition to be in contact
     * @param ActiveCheckFactor: The proportion of the length of the geometry that is going to be taking into account to chechk if the node is active or inactive
     * @return ConditionPointers: The vector containing the pointer to the conditions of interest
     */
    
    void MortarContainerFiller(
        Condition::Pointer & pCond_1,
        const Condition::Pointer & pCond_2,
        std::vector<contact_container> *& ConditionPointers,
        const double ActiveCheckFactor
        )
    {
        SearchUtilities::ContactContainerFiller(ConditionPointers, pCond_1, pCond_2, 
                        pCond_1->GetValue(NORMAL), pCond_2->GetValue(NORMAL), ActiveCheckFactor); 
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     */
    
    void ComputeNodesMeanNormal()
    {
        ContactUtilities::ComputeNodesMeanNormalModelPart(mrDestinationModelPart);
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
        return "DeprecatedTreeContactSearch";
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

    void ResetContactOperators(ModelPart & rModelPart)
    {
        ConditionsArrayType& pCond = rModelPart.Conditions();
    
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            if (itCond->Is(ACTIVE) == true)
            {
                itCond->Set(ACTIVE, false);
                
                std::vector<contact_container> * ConditionPointers = itCond->GetValue(CONTACT_CONTAINERS);
                
                if (ConditionPointers != NULL)
                {
                    for (unsigned int i = 0; i < ConditionPointers->size();i++)
                    {
                        (*ConditionPointers)[i].clear();
                    } 
                    
                    ConditionPointers->clear();
//                     ConditionPointers->reserve(mallocation); 
                }
//                 delete ConditionPointers;
//                 itCond->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
            }
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
  
    ModelPart& mrOriginModelPart;             // The original model part
    ModelPart& mrDestinationModelPart;        // The destination model part
    unsigned int mBucketSize;                 // Bucket size for kd-tree
    PointVector mPointListDestination;        // A list that contents the all the points (from nodes) from the modelpart 
    const unsigned int mdimension;            // Dimension size of the space
    const unsigned int mallocation;           // Allocation size for the vectors and max number of potential results
    
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

}; // Class DeprecatedTreeContactSearch

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
                                  DeprecatedTreeContactSearch& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DeprecatedTreeContactSearch& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
