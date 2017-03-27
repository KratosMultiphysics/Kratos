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
#include <iostream>
#include <vector>
#include "boost/smart_ptr.hpp"

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product
#include "custom_utilities/contact_utilities.h"

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

#if !defined(POINT_ITEM_DEFINED )
#define  POINT_ITEM_DEFINED
/** @brief Custom Point container to be used by the mapper
 */
class PointItem: public Point<3>
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of PointItem
    KRATOS_CLASS_POINTER_DEFINITION( PointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointItem():
        Point<3>()
    {
    }

    PointItem(array_1d<double, 3> Coords):
        Point<3>(Coords)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Condition::Pointer Cond,
        double Radius
    ):
        Point<3>(Coords),
        mpOriginCond(Cond),
        mRadius(Radius)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Node<3>::Pointer Node
    ):
        Point<3>(Coords),
        mpOriginNode(Node)
    {}

    ///Copy constructor  (not really required)
    PointItem(const PointItem& rhs):
        Point<3>(rhs),
        mpOriginCond(rhs.mpOriginCond),
        mpOriginNode(rhs.mpOriginNode),
        mRadius(rhs.mRadius)
    {
    }

    /// Destructor.
    // ~PointItem();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Returns the point
     * @return The point
     */
    Point<3> GetPoint()
    {
        Point<3> Point(this->Coordinates());
        
        return Point;
    }
    
    /**
     * Set the point
     * @param The point
     */
    void SetPoint(Point<3> Point)
    {
        this->Coordinates() = Point.Coordinates();
    }
    
    /**
     * Returns the radius of the condition
     * @return mRadius: The radius of the condition
     */
    double GetRadius()
    {
        return mRadius;
    }
    
    /**
     * Sets the radius of the condition
     * @param Radius: The radius of the condition
     */
    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    /**
     * Sets the condition associated to the point
     * @param Cond: The pointer to the condition
     */

    void SetCondition(Condition::Pointer Cond)
    {
        mpOriginCond = Cond;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginCond: The pointer to the condition associated to the point
     */

    Condition::Pointer GetCondition()
    {
        return mpOriginCond;
    }
    
    /**
     * Sets the node associated to the point
     * @param Node: The pointer to the node associated to the point
     */

    void SetNode(Node<3>::Pointer Node)
    {
        mpOriginNode = Node;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginNode: The pointer to the node associated to the point
     */

    Node<3>::Pointer GetNode()
    {
        return mpOriginNode;
    }

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

    Condition::Pointer mpOriginCond; // Condition pointer
    Node<3>::Pointer   mpOriginNode; // Node pointer
    double                  mRadius; // Radius         
    array_1d<double, 3>     mNormal; // Normal vector      

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PointItem    
#endif // POINT_ITEM_DEFINED  defined 
    
/** \brief TreeContactSearch
 * This utilitiy has as objective to create the contact conditions. The conditions that can be created are:
 * * Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments (En principio)
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree (En principio no)
 * To consider autocontact use the same model_part as origin and destination (En principio, preguntar a Riccardo)
 */

class TreeContactSearch
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

    /// Pointer definition of TreeContactSearch
    // KRATOS_CLASS_POINTER_DEFINITION( TreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{
    
    // Class Constructor
    // WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
    // Use an InterfacePreprocess object to create such a model part from a regular one:
    // InterfaceMapper = InterfacePreprocess()
    // InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
    TreeContactSearch(
        ModelPart & rOriginModelPart,
        ModelPart & rDestinationModelPart,
        const unsigned int allocation_size
        )
    :mrOriginModelPart(rOriginModelPart),
     mrDestinationModelPart(rDestinationModelPart),
     mBucketSize(4),
     mdimension(rOriginModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension()),
     mallocation(allocation_size)
    {
//         // Destination model part
//         ModelPartSetter(mrDestinationModelPart, false, true, false);
        
        NodesArrayType& pNode = mrDestinationModelPart.Nodes();
            
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->Set(SLAVE, true);  
        }

//         // Origin model part
//         ModelPartSetter(mrOriginModelPart, false, false, true); 
        
        ConditionsArrayType& pCond = mrOriginModelPart.Conditions();
            
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            itCond->Set(MASTER, true);  
        }
    }
        
    void ModelPartSetter(
        ModelPart & rModelPart,
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
    
    void InitializeALMFrictionlessMortarConditions(
        const double rActiveCheckFactor,
        const int rIntegrationOrder
        )
    {
        // Destination model part
        InitializeALMFrictionlessConditions(mrDestinationModelPart, rActiveCheckFactor, rIntegrationOrder);
        
        // Origin model part
        InitializeALMFrictionlessConditions(mrOriginModelPart, rActiveCheckFactor, rIntegrationOrder);
    }
    
    /**
     * This function initializes nodes 
     */
    
    void InitializeNodes(ModelPart & rModelPart);
    
    /**
     * This function initializes ALM frictionless conditions
     */
        
    void InitializeALMFrictionlessConditions(
        ModelPart & rModelPart, 
        const double rActiveCheckFactor,
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
            if (cond_it->GetProperties().Has(INTEGRATION_ORDER_CONTACT) == false)
            {
                cond_it->GetProperties().SetValue(INTEGRATION_ORDER_CONTACT, rIntegrationOrder);
            }
        }
    }
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearALMFrictionlessMortarConditions()
    {
        // Destination model part
        TotalClearALMFrictionlessConditions(mrDestinationModelPart);   
    }
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void PartialClearALMFrictionlessMortarConditions()
    {
        // Destination model part
        PartialClearALMFrictionlessConditions(mrDestinationModelPart); 
    }
    
    /**
     * This function clears conditions already created 
     */
    
    void TotalClearALMFrictionlessConditions(ModelPart & rModelPart)
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
                itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
            }
        }  
    }
    
    /**
     * This function clears partially the conditions already created 
     */
    
    void PartialClearALMFrictionlessConditions(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            if (itNode->Is(ACTIVE) == false)
            {
                itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
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
     * This function creates a node list
     */
    
    void CreatePointListNodes(
        ModelPart & rModelPart, 
        PointVector & PoinList
        )
    {
        array_1d<double, 3> Coord = ZeroVector(3); // Will store the coordinates 

        NodesArrayType& pNode               = rModelPart.Nodes();
        NodesArrayType::iterator node_begin = pNode.ptr_begin();
        NodesArrayType::iterator node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = node_begin; node_it!=node_end; node_it++)
        {
            noalias(Coord) = node_it->Coordinates();
            
            PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); 
            (PoinList).push_back(pP);
        }
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
            Point<3> Center;
            double Radius = ContactUtilities::CenterAndRadius(pCond, Center); 
            PointItem::Pointer pPoint = PointItem::Pointer(new PointItem(Center, pCond, Radius));
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
     * This function updates a node list
     */
    
    void UpdatePointListNodes(
        ModelPart & rModelPart, 
        PointVector & PoinList
        )
    {
        // TODO: Add this in the future
    }

    /**
     * This function updates a condition list
     */
    
    void UpdatePointListConditions(
        ModelPart & rModelPart, 
        PointVector & PoinList
        )
    {
        // TODO: Think how to parallel this!!!!
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_end   = pCond.ptr_end();
        
        unsigned int index = 0;
        for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
        {
            const Condition::Pointer pCond = (*cond_it.base());
            Point<3> Center;
            const double Radius = ContactUtilities::CenterAndRadius(pCond, Center); 
            PointItem::Pointer & pPoint = PoinList[index];
            pPoint->SetCondition(pCond);
            pPoint->SetRadius(Radius);
            pPoint->SetPoint(Center);
            index += 1;
        }
    }
    
    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
    
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

                NumberPointsFound = Tree_points.SearchInRadius(Center, SearchRadius, PointsFound.begin(), PointsDistances.begin(), mallocation);
            }
//             else if (type_search == 1) // TODO: Complete search in bounding box
//             {
//                 Point<3> MinPoint, MaxPoint;
//                 CondOri->GetGeometry().BoundingBox(MinPoint, MaxPoint);
//                 NumberPointsFound= Tree_conds.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), PointsDistances.begin(), mallocation);
//             }
//             else if (type_search == 1) // TODO: Complete search in k-DOP
//             {
//             }
            else
            {
                KRATOS_ERROR << " The type search declared does not exist!!!!. type_search = " << type_search << std::endl;
            }
            
            if (NumberPointsFound > 0)
            {
    //             KRATOS_WATCH(NumberPointsFound); 
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
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
        

    void CreateALMFrictionlessMortarConditions(
        const double SearchFactor,
        const int type_search
    )
    {
        TotalClearALMFrictionlessMortarConditions(); // Clear the conditions
        UpdateMortarConditions(SearchFactor, type_search);
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
                std::cout << "Node: " << node_it->Id() << " is active. SLAVE: " << node_it->Is(SLAVE) << std::endl;
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
        
//         #pragma omp parallel for // NOTE: Be careful, if you change something with get value over the nodes iteraring in the conditions this will not work in OpenMP 
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
                            const bool cond_active = ContactUtilities::ContactContainerFiller(itCond->GetGeometry(),    (aux_contact_container.condition)->GetGeometry(), 
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
        ContactUtilities::ContactContainerFiller(ConditionPointers, pCond_1, pCond_2, 
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
