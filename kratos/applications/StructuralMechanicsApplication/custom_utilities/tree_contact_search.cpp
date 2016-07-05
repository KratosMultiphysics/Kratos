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
    AuxConstructor(mrDestinationModelPart, true);
    
    // Origin model part
    AuxConstructor(mrOriginModelPart, false);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::AuxConstructor(
    ModelPart & rModelPart,
    const bool rSlaveMaster
    ) 
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->Set( ACTIVE, false ); // NOTE: It is supposed to be already false, just in case   
        cond_it->Set( SLAVE,   rSlaveMaster);
        cond_it->Set( MASTER, !rSlaveMaster);
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

void TreeContactSearch::InitializeMortarConditions()
{
    // Destination model part
    InitializeConditions(mrDestinationModelPart);
    
    // Origin model part
    InitializeConditions(mrOriginModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeNodes(ModelPart & rModelPart)
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeConditions(ModelPart & rModelPart)
{
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
//         cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
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
    
    // Origin model part
    ClearConditions(mrOriginModelPart);
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
        if (cond_it->Is(ACTIVE))
        {
            cond_it->Set(ACTIVE, false);
            
            std::vector<contact_container> *& ConditionPointers = cond_it->GetValue(CONTACT_CONTAINERS);
            
            for (unsigned int i =0; i< ConditionPointers->size();i++)
            {
//                 delete(&((*ConditionPointers)[i].condition));
                (*ConditionPointers)[i].clear();
            } 
              
            ConditionPointers->clear();
//             ConditionPointers->reserve(mallocation); 
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTN()
{
    // Destination model part
    CreatePointListNodes(mrDestinationModelPart, mPointListDestination);
    
    // Origin model part
    CreatePointListNodes(mrOriginModelPart, mPointListOrigin);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTS()
{
    // Destination model part
    CreatePointListNodes(mrDestinationModelPart, mPointListDestination);
    
    // Origin model part
    CreatePointListConditions(mrOriginModelPart, mPointListOrigin);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListMortar()
{
    // Destination model part
    CreatePointListConditions(mrDestinationModelPart, mPointListDestination);
    
    // Origin model part
    CreatePointListConditions(mrOriginModelPart, mPointListOrigin);
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
        const Condition::Pointer Cond = *cond_it.base();
        Point<3> Center;
        double Radius;
        ContactUtilities::CenterAndRadius(Cond, Center, Radius, mdimension); 
        PointItem::Pointer pPoint = PointItem::Pointer(new PointItem(Center, *(cond_it.base()), Radius));
        (PoinList).push_back(pPoint);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdatePointListMortar()
{
    // Destination model part
    UpdatePointListConditions(mrDestinationModelPart, mPointListDestination);
    
    // Origin model part
    UpdatePointListConditions(mrOriginModelPart, mPointListOrigin);
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
        const Condition::Pointer Cond = *cond_it.base();
        Point<3> Center;
        double Radius;
        ContactUtilities::CenterAndRadius(Cond, Center, Radius, mdimension); 
        PointItem::Pointer & pPoint = PoinList[index];
        pPoint->SetRadius(Radius);
        pPoint->SetPoint(Center);
        index += 1;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreateNTNConditions(
    const double SearchFactor,
    const unsigned int MaxNumberResults,
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
    const unsigned int MaxNumberResults,
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
    const unsigned int MaxNumberResults,
    const int type_search,
    const int integration_order
) 
{
    // Initialize values
    PointVector PointsFound(MaxNumberResults);
    std::vector<double> PointsDistances(MaxNumberResults);
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
    
    for(unsigned int cond_it = 0; cond_it < mPointListOrigin.size(); cond_it++)
    {
        if (type_search == 0)
        {
            double SearchRadius = SearchFactor * mPointListOrigin[cond_it]->GetRadius();
            Point<3> Center = mPointListOrigin[cond_it]->GetPoint();
            NumberPointsFound = Tree_points.SearchInRadius(Center, SearchRadius, PointsFound.begin(), PointsDistances.begin(), MaxNumberResults);
        }
//         else if (type_search == 1) // TODO: Complete search in bounding box
//         {
//             Point<3> MinPoint, MaxPoint;
//             CondOri->GetGeometry().BoundingBox(MinPoint, MaxPoint);
//             NumberPointsFound= Tree_conds.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), PointsDistances.begin(), MaxNumberResults);
//         }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error, " The type search declared does not exist!!!!. type_search = ", type_search );
        }
        
        if (NumberPointsFound > 0)
        {
            Condition::Pointer pCondOrigin = mPointListOrigin[cond_it]->GetCondition();
            
//             KRATOS_WATCH(NumberPointsFound); 
            for(unsigned int i = 0; i < NumberPointsFound; i++)
            {   
                Condition::Pointer pCondDestination = PointsFound[i]->GetCondition();
                std::vector<contact_container> *& ConditionPointersDestination = pCondDestination->GetValue(CONTACT_CONTAINERS);
                int & DestCondIntegrationOrder = pCondDestination->GetValue(INTEGRATION_ORDER_CONTACT);
                DestCondIntegrationOrder = integration_order;
                
                // Set the corresponding flags
                pCondDestination->Set(ACTIVE, true);

                MortarContainerFiller(mPointListOrigin[cond_it], PointsFound[i], pCondDestination, pCondOrigin, ConditionPointersDestination, IntegrationOrder, false);
            }
        }
    }
    
    // FIXME
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
    
    ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        if (cond_it->Is(ACTIVE))
        {
            std::vector<contact_container> *& ConditionPointersOrigin = cond_it->GetValue(CONTACT_CONTAINERS);
            KRATOS_WATCH(ConditionPointersOrigin->size());
            
            for (unsigned int i = 0; i < ConditionPointersOrigin->size(); i++)
            {
                (*ConditionPointersOrigin)[i].print();
            } 
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::MortarContainerFiller(
        const PointType::Pointer PointOfList,
        const PointType::Pointer PointFound,
        const Condition::Pointer & pCond_1,
        const Condition::Pointer & pCond_2,
        std::vector<contact_container> *& ConditionPointers,
        const IntegrationMethod & IntegrationOrder,
        const bool orig_dest
        )
{
    contact_container contact_container;
    
    contact_container.condition = &* pCond_2;
    
    // Define the normal to the contact
    Point<3> ContactPoint;
    if (orig_dest == true)
    {
        ContactPoint = PointFound->GetPoint();
    }
    else
    {
        ContactPoint = PointOfList->GetPoint();
    }
    
    ContactUtilities::ContactContainerFiller(contact_container, ContactPoint, pCond_1->GetGeometry(), pCond_2->GetGeometry(), 
                                             pCond_1->GetValue(NORMAL), pCond_2->GetValue(NORMAL), IntegrationOrder);
        
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
