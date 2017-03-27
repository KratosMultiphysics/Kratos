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

// System includes

// External includes

// System includes
# include "deprecated_tree_contact_search.h"

// TODO: Not clear everything from one step to the other, just check what can you add (this way the process is simplified)

namespace Kratos
{
/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/    

// Class Constructor
// WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
// Use an InterfacePreprocess object to create such a model part from a regular one:
// InterfaceMapper = InterfacePreprocess()
// InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
DeprecatedTreeContactSearch::DeprecatedTreeContactSearch(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::ModelPartSetter(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::InitializeMortarConditions(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::InitializeMeshTyingMortarConditions(
    const double rActiveCheckFactor,
    const int rIntegrationOrder
    )
{
    // Destination model part
    InitializeMeshTyingConditions(mrDestinationModelPart, rActiveCheckFactor, rIntegrationOrder);
    
    // Origin model part
    InitializeMeshTyingConditions(mrOriginModelPart, rActiveCheckFactor, rIntegrationOrder);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::InitializeNodes(ModelPart & rModelPart)
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::InitializeConditions(
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
//     
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
//         cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
        cond_it->GetProperties().SetValue(ACTIVE_CHECK_FACTOR, rActiveCheckFactor);
        cond_it->GetProperties().SetValue(NORMAL_AUGMENTATION_FACTOR,  rAugmentationNormal);
        cond_it->GetProperties().SetValue(TANGENT_AUGMENTATION_FACTOR, rAugmentationTangent);
        if (cond_it->GetProperties().Has(INTEGRATION_ORDER_CONTACT) == false)
        {
            cond_it->GetProperties().SetValue(INTEGRATION_ORDER_CONTACT, rIntegrationOrder);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::InitializeMeshTyingConditions(
    ModelPart & rModelPart,
    const double rActiveCheckFactor,
    const int rIntegrationOrder
    )
{   
    ConditionsArrayType& pCond  = rModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCond.ptr_end();
//     
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
//         cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
        cond_it->GetProperties().SetValue(ACTIVE_CHECK_FACTOR, rActiveCheckFactor);
        if (cond_it->GetProperties().Has(INTEGRATION_ORDER_CONTACT) == false)
        {
            cond_it->GetProperties().SetValue(INTEGRATION_ORDER_CONTACT, rIntegrationOrder);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearMortarConditions()
{
    // Destination model part
    TotalClearConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearMeshTyingMortarScalarConditions()
{
    // Destination model part
    TotalClearMeshTyingScalarConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearMeshTyingMortarComponentsConditions()
{
    // Destination model part
    TotalClearMeshTyingComponentsConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearMortarConditions()
{
    // Destination model part
    PartialClearConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearMeshTyingMortarScalarConditions()
{
    // Destination model part
    PartialClearMeshTyingScalarConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearMeshTyingMortarComponentsConditions()
{
    // Destination model part
    PartialClearMeshTyingComponentsConditions(mrDestinationModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearConditions(ModelPart & rModelPart)
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearMeshTyingScalarConditions(ModelPart & rModelPart)
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
            itNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER, 0) = 0.0;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::TotalClearMeshTyingComponentsConditions(ModelPart & rModelPart)
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearConditions(ModelPart & rModelPart)
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearMeshTyingScalarConditions(ModelPart & rModelPart)
{
    NodesArrayType& pNode = rModelPart.Nodes();
    
    auto numNodes = pNode.end() - pNode.begin();
    
    #pragma omp parallel for 
    for(unsigned int i = 0; i < numNodes; i++) 
    {
        auto itNode = pNode.begin() + i;
        if (itNode->Is(ACTIVE) == false)
        {
            itNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER, 0) = 0.0;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::PartialClearMeshTyingComponentsConditions(ModelPart & rModelPart)
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::CreatePointListMortar()
{
    // Destination model part
    CreatePointListConditions(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::CreatePointListNodes(
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

void DeprecatedTreeContactSearch::CreatePointListConditions(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::UpdatePointListMortar()
{
    // Destination model part
    UpdatePointListConditions(mrDestinationModelPart, mPointListDestination);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::UpdatePointListNodes(
    ModelPart & rModelPart, 
    PointVector & PoinList
    )
{
    // TODO: Add this in the future
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::UpdatePointListConditions(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::CreateMortarConditions(
    const double SearchFactor,
    const int type_search
) 
{
    TotalClearMortarConditions(); // Clear the conditions
    UpdateMortarConditions(SearchFactor, type_search);
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::UpdateMortarConditions( // TODO: Change everything, using the slave as reference isntead of the master
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
//         else if (type_search == 1) // TODO: Complete search in bounding box
//         {
//             Point<3> MinPoint, MaxPoint;
//             CondOri->GetGeometry().BoundingBox(MinPoint, MaxPoint);
//             NumberPointsFound= Tree_conds.SearchInBox(MinPoint, MaxPoint, PointsFound.begin(), PointsDistances.begin(), mallocation);
//         }
//         else if (type_search == 1) // TODO: Complete search in k-DOP
//         {
//         }
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::CheckMortarConditions()
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::ClearAllInactivePairs(ModelPart & rModelPart)
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

/***********************************************************************************/
/***********************************************************************************/

bool DeprecatedTreeContactSearch::CheckCondition(
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

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::MortarContainerFiller(
        Condition::Pointer & pCondDestination,
        const Condition::Pointer & pCondOrigin,
        std::vector<contact_container> *& ConditionPointers,
        const double ActiveCheckFactor
        )
{
    ContactUtilities::ContactContainerFiller(ConditionPointers, pCondDestination, pCondOrigin, 
                      pCondDestination->GetValue(NORMAL), pCondOrigin->GetValue(NORMAL), ActiveCheckFactor); 
}

/***********************************************************************************/
/***********************************************************************************/

void DeprecatedTreeContactSearch::ComputeNodesMeanNormal()
{
    ContactUtilities::ComputeNodesMeanNormalModelPart(mrDestinationModelPart);
}

/************************************** ACCESS *************************************/
/***********************************************************************************/
    
}  // namespace Kratos.
