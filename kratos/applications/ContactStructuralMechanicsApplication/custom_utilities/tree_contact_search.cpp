// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

// System includes

// External includes

// System includes
# include "tree_contact_search.h"

namespace Kratos
{
/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/    

// Class Constructor
// WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
// Use an InterfacePreprocess object to create such a model part from a regular one:
// InterfaceMapper = InterfacePreprocess()
// InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
TreeContactSearch::TreeContactSearch(
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
    // Destination model part
    AuxConstructor(mrDestinationModelPart, false, false);
    
    // Origin model part
    AuxConstructor(mrOriginModelPart, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::AuxConstructor(
    ModelPart & rModelPart,
    const bool rActive,
//     const bool rSlave,
    const bool rMaster
    ) 
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->Set( ACTIVE, rActive ); // NOTE: It is supposed to be already false, just in case   
//         cond_it->Set( SLAVE,  rSlave);
        cond_it->Set( MASTER, rMaster);
    }
    
    NodesArrayType& pNode               = rModelPart.Nodes();
    NodesArrayType::iterator node_begin = pNode.ptr_begin();
    NodesArrayType::iterator node_end   = pNode.ptr_end();
    
    for(NodesArrayType::iterator node_it = node_begin; node_it!=node_end; node_it++)
    {
        node_it->Set( ACTIVE, rActive );  // NOTE: It is supposed to be already false, just in case   
    }
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

TreeContactSearch::~TreeContactSearch() {}


/************************************* OPERATIONS **********************************/
/***********************************************************************************/

void TreeContactSearch::InitializeNTNConditions()
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeNTSConditions()
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeMortarConditions(
    const double rActiveCheckFactor,
    const double rConstantActInact
    )
{
    // Destination model part
    InitializeConditions(mrDestinationModelPart, rActiveCheckFactor, rConstantActInact);
    
    // Origin model part
    InitializeConditions(mrOriginModelPart, rActiveCheckFactor, rConstantActInact);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeNodes(ModelPart & rModelPart)
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeConditions(
    ModelPart & rModelPart,
    const double rActiveCheckFactor,
    const double rConstantActInact
    )
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
//         cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
        cond_it->GetProperties().SetValue(ACTIVE_CHECK_FACTOR, rActiveCheckFactor);
        cond_it->GetProperties().SetValue(CONSTANT_ACT_INACT,  rConstantActInact);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ClearNTNConditions()
{    
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ClearNTSConditions()
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ClearMortarConditions()
{
    // Destination model part
    ClearConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ClearConditions(ModelPart & rModelPart)
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        if (cond_it->Is(ACTIVE) == true)
        {
            cond_it->Set(ACTIVE, false);
            
            std::vector<contact_container> * ConditionPointers = cond_it->GetValue(CONTACT_CONTAINERS);
//             std::vector<contact_container> *& ConditionPointers = cond_it->GetValue(CONTACT_CONTAINERS);
            
            if (ConditionPointers != NULL)
            {
                for (unsigned int i = 0; i < ConditionPointers->size();i++)
                {
                    (*ConditionPointers)[i].clear();
                } 
                
                ConditionPointers->clear();
                
    //             ConditionPointers->reserve(mallocation); 
            }
        
//             delete ConditionPointers;
//             cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
        }
    }
    
    NodesArrayType& pNode               = rModelPart.Nodes();
    NodesArrayType::iterator node_begin = pNode.ptr_begin();
    NodesArrayType::iterator node_end   = pNode.ptr_end();
    
    for(NodesArrayType::iterator node_it = node_begin; node_it!=node_end; node_it++)
    {
        if (node_it->Is(ACTIVE) == true)
        {
            node_it->Set( ACTIVE, false );
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTN()
{
    // Destination model part
    CreatePointListNodes(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTS()
{
    // Destination model part
    CreatePointListNodes(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListMortar()
{
    // Destination model part
    CreatePointListConditions(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNodes(
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

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListConditions(
    ModelPart & rModelPart, 
    PointVector & PoinList
    )
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        const Condition::Pointer pCond = (*cond_it.base());
        Point<3> Center;
        double Radius;
        ContactUtilities::CenterAndRadius(pCond, Center, Radius, mdimension); 
        PointItem::Pointer pPoint = PointItem::Pointer(new PointItem(Center, pCond, Radius));
        (PoinList).push_back(pPoint);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdatePointListMortar()
{
    // Destination model part
    UpdatePointListConditions(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdatePointListNodes(
    ModelPart & rModelPart, 
    PointVector & PoinList
    )
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdatePointListConditions(
    ModelPart & rModelPart, 
    PointVector & PoinList
    )
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    unsigned int index = 0;
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        const Condition::Pointer pCond = (*cond_it.base());
        Point<3> Center;
        double Radius;
        ContactUtilities::CenterAndRadius(pCond, Center, Radius, mdimension); 
        PointItem::Pointer & pPoint = PoinList[index];
        pPoint->SetCondition(pCond);
        pPoint->SetRadius(Radius);
        pPoint->SetPoint(Center);
        index += 1;
     }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreateNTNConditions(
    const double SearchFactor,
    const int type_search,
    const int integration_order
) 
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreateNTSConditions(
    const double SearchFactor,
    const int type_search,
    const int integration_order
) 
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreateMortarConditions(
    const double SearchFactor,
    const int type_search,
    const int integration_order
) 
{
    // Initialize values
    PointVector PointsFound(mallocation);
    std::vector<double> PointsDistances(mallocation);
    unsigned int NumberPointsFound = 0;
    ClearMortarConditions(); // Clear the conditions
    
    IntegrationMethod IntegrationOrder;

    if (integration_order == 1)
    {
        IntegrationOrder = GeometryData::GI_GAUSS_1;
    }
    else if (integration_order == 2)
    {
        IntegrationOrder = GeometryData::GI_GAUSS_2;
    }
    else if (integration_order == 3)
    {
        IntegrationOrder = GeometryData::GI_GAUSS_3;
    }
    else if (integration_order == 4)
    {
        IntegrationOrder = GeometryData::GI_GAUSS_4;
    }
    else if (integration_order == 5)
    {
        IntegrationOrder = GeometryData::GI_GAUSS_5;
    }
    else
    {
        std::cout << "The number of integration points is not defined.  integration_order: "<< integration_order << std::endl;
        std::cout << "Taking default number of integration points" << std::endl;
        IntegrationOrder = mrOriginModelPart.ConditionsBegin()->GetGeometry().GetDefaultIntegrationMethod();
    }
    
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
            double SearchRadius;
            ContactUtilities::CenterAndRadius(pCondOrigin, Center, SearchRadius, mdimension); 
            SearchRadius *= SearchFactor;

            NumberPointsFound = Tree_points.SearchInRadius(Center, SearchRadius, PointsFound.begin(), PointsDistances.begin(), mallocation);
        }
//         else if (type_search == 1) // TODO: Complete search in bounding box
//         {
//             Point<3> MinPoint, MaxPoint;
//             CondOri->GetGeometry().BoundingBox(MinPoint, MaxPoint);
//             NumberPointsFound= Tree_conds.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), PointsDistances.begin(), mallocation);
//         }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error, " The type search declared does not exist!!!!. type_search = ", type_search );
        }
        
        if (NumberPointsFound > 0)
        {
//             KRATOS_WATCH(NumberPointsFound); 
            for(unsigned int i = 0; i < NumberPointsFound; i++)
            {   
                Condition::Pointer pCondDestination = PointsFound[i]->GetCondition();
                
                if (pCondDestination != pCondOrigin) // Avoiding "auto self-contact"
                {
                    std::vector<contact_container> *& ConditionPointersDestination = pCondDestination->GetValue(CONTACT_CONTAINERS);
                    int & DestCondIntegrationOrder = pCondDestination->GetValue(INTEGRATION_ORDER_CONTACT);
                    DestCondIntegrationOrder = integration_order;
                    
                    // Set the corresponding flags
                    pCondDestination->Set(ACTIVE, true); 
                    const double ActiveCheckFactor = pCondDestination->GetProperties().GetValue(ACTIVE_CHECK_FACTOR);

                    MortarContainerFiller(pCondOrigin->GetGeometry().Center(), PointsFound[i], pCondDestination, pCondOrigin, ConditionPointersDestination, ActiveCheckFactor, IntegrationOrder, false);
                }
            }
        }
    }
    
    // Calculate the mean of the normal in all the nodes
    ComputeNodesMeanNormal();
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CheckMortarConditions()
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

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::MortarContainerFiller(
        const Point<3>& OriginPoint,
        const PointType::Pointer PointFound,
        const Condition::Pointer & pCond_1,
        const Condition::Pointer & pCond_2,
        std::vector<contact_container> *& ConditionPointers,
        const double ActiveCheckFactor,
        const IntegrationMethod & IntegrationOrder,
        const bool orig_dest
        )
{
    contact_container contact_container;
    
    contact_container.condition = pCond_2;
    
    // Define the normal to the contact
    Point<3> ContactPoint;
    if (orig_dest == true)
    {
        ContactPoint = PointFound->GetPoint();
    }
    else
    {
        ContactPoint = OriginPoint;
    }
    
    ContactUtilities::ContactContainerFiller(contact_container, ContactPoint, pCond_1->GetGeometry(), pCond_2->GetGeometry(), 
                                             pCond_1->GetValue(NORMAL), pCond_2->GetValue(NORMAL), ActiveCheckFactor, IntegrationOrder); 
        
    ConditionPointers->push_back(contact_container);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ComputeNodesMeanNormal()
{
    ContactUtilities::ComputeNodesMeanNormalModelPart(mrDestinationModelPart);
}

/************************************** ACCESS *************************************/
/***********************************************************************************/
    
}  // namespace Kratos.
