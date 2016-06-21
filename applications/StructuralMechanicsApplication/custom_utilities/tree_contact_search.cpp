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
    ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCondDestination.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCondDestination.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->Set( ACTIVE, false ); // NOTE: It is supposed to be already false, just in case
        cond_it->Set( MASTER, false ); // NOTE: It is supposed to be already false, just in case
        cond_it->Set( SLAVE,  true  );      
    }
    
    ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->Set( ACTIVE, false ); // NOTE: It is supposed to be already false, just in case
        cond_it->Set( MASTER, true  );
        cond_it->Set( SLAVE,  false ); // NOTE: It is supposed to be already false, just in case   
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
    ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCondDestination.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCondDestination.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
        cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
    }
    
    ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();   
        cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation);    
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
    ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCondDestination.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCondDestination.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        if (cond_it->Is(ACTIVE))
        {
            cond_it->Set(ACTIVE, false);
            
            std::vector<contact_container> *& ConditionPointersDestination = cond_it->GetValue(CONTACT_CONTAINERS);
            
            for (unsigned int i =0; i< ConditionPointersDestination->size();i++)
            {
                delete(&((*ConditionPointersDestination)[i].condition));
            } 
            
            ConditionPointersDestination->clear();
            
            delete cond_it->GetValue(CONTACT_CONTAINERS); 
            cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>(); 
            cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
        }
    }
    
    ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        if (cond_it->Is(ACTIVE))
        {
            cond_it->Set(ACTIVE, false); 
            
            std::vector<contact_container> *& ConditionPointersOrigin = cond_it->GetValue(CONTACT_CONTAINERS);
            
            for (unsigned int i =0; i< ConditionPointersOrigin->size();i++)
            {
                delete(&((*ConditionPointersOrigin)[i].condition));
            } 
            
            ConditionPointersOrigin->clear();
            
            delete cond_it->GetValue(CONTACT_CONTAINERS); 
            cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
            cond_it->GetValue(CONTACT_CONTAINERS)->reserve(mallocation); 
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTN()
{
    array_1d<double, 3> Coord; // Will store the coordinates 
    Coord = ZeroVector(3);
    
    NodesArrayType& pNodeDestination  = mrDestinationModelPart.Nodes();
    NodesArrayType::iterator it_begin = pNodeDestination.ptr_begin();
    NodesArrayType::iterator it_end   = pNodeDestination.ptr_end();
    
    for(NodesArrayType::iterator node_it = it_begin; node_it!=it_end; node_it++)
    {
        Coord[0] = node_it->X();
        Coord[1] = node_it->Y();
        Coord[2] = node_it->Z();

        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); 
        (mPointListDestination).push_back(pP);
    }
    
    NodesArrayType& pNodeOrigin  = mrOriginModelPart.Nodes();
    it_begin = pNodeOrigin.ptr_begin();
    it_end   = pNodeOrigin.ptr_end();
    
    for(NodesArrayType::iterator node_it = it_begin; node_it!=it_end; node_it++)
    {
        Coord[0] = node_it->X();
        Coord[1] = node_it->Y();
        Coord[2] = node_it->Z();
        
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); 
        (mPointListOrigin).push_back(pP);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListNTS()
{
    array_1d<double, 3> Coord; // Will store the coordinates 
    Coord = ZeroVector(3);
    
    NodesArrayType& pNode               = mrDestinationModelPart.Nodes();
    NodesArrayType::iterator node_begin = pNode.ptr_begin();
    NodesArrayType::iterator node_end   = pNode.ptr_end();
    
    for(NodesArrayType::iterator node_it = node_begin; node_it!=node_end; node_it++)
    {
        Coord[0] = node_it->X();
        Coord[1] = node_it->Y();
        Coord[2] = node_it->Z();
        
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); 
        (mPointListDestination).push_back(pP);
    }
    
    ConditionsArrayType& pCond               = mrOriginModelPart.Conditions();
    ConditionsArrayType::iterator cond_begin = pCond.ptr_begin();
    ConditionsArrayType::iterator cond_end   = pCond.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = cond_begin; cond_it!=cond_end; cond_it++)
    {
        const Condition::Pointer CondOri = *cond_it.base();
        Point<3> Center;
        double Radius;
        array_1d<double, 3> Normal;
        StructuralMechanicsMathUtilities::CenterAndRadius(CondOri, Center, Radius, Normal, mdimension); 
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Center, *(cond_it.base()), Radius, Normal));
        (mPointListOrigin).push_back(pP);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListMortar()
{
    ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
    ConditionsArrayType::iterator it_begin = pCondDestination.ptr_begin();
    ConditionsArrayType::iterator it_end   = pCondDestination.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        const Condition::Pointer CondOri = *cond_it.base();
        Point<3> Center;
        double Radius;
        array_1d<double, 3> Normal;
        StructuralMechanicsMathUtilities::CenterAndRadius(CondOri, Center, Radius, Normal, mdimension); 
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Center, *(cond_it.base()), Radius, Normal));
        (mPointListDestination).push_back(pP);
    }
    
    ConditionsArrayType& pCondOrigin  = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        const Condition::Pointer CondOri = *cond_it.base();
        Point<3> Center;
        double Radius;
        array_1d<double, 3> Normal;
        StructuralMechanicsMathUtilities::CenterAndRadius(CondOri, Center, Radius, Normal, mdimension); 
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Center, *(cond_it.base()), Radius, Normal));
        (mPointListOrigin).push_back(pP);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreateNTNConditions(
    const double SearchFactor,
    const unsigned int MaxNumberResults,
    const int type_search,
    const bool bidirectional,
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
    const bool bidirectional,
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
    const bool bidirectional,
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

            if (bidirectional == true)
            {
                // Set the corresponding flags
                pCondOrigin->Set(ACTIVE, true);
            }

            std::vector<contact_container> *& ConditionPointersOrigin = pCondOrigin->GetValue(CONTACT_CONTAINERS);
            pCondOrigin->GetValue(INTEGRATION_ORDER_CONTACT) = IntegrationOrder;
            
//             KRATOS_WATCH(NumberPointsFound); 
            for(unsigned int i = 0; i < NumberPointsFound; i++)
            {   
                Condition::Pointer pCondDestination = PointsFound[i]->GetCondition();
                std::vector<contact_container> *& ConditionPointersDestination = pCondDestination->GetValue(CONTACT_CONTAINERS);
                pCondDestination->GetValue(INTEGRATION_ORDER_CONTACT) = IntegrationOrder;
                
                if (bidirectional == true)
                {
                    MortarContainerFiller(mPointListOrigin[cond_it], PointsFound[i], pCondOrigin, pCondDestination, ConditionPointersOrigin, IntegrationOrder, true);
                }
                
                // Set the corresponding flags
                pCondDestination->Set(ACTIVE, true);

                MortarContainerFiller(mPointListOrigin[cond_it], PointsFound[i], pCondDestination, pCondOrigin, ConditionPointersDestination, IntegrationOrder, false);
            }
        }
    }
    
    // FIXME
    // Calculate the mean of the normal in all the nodes
    ComputeNodesMeanNormal(bidirectional);
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
    array_1d<double, 3> & contact_normal = pCond_1->GetValue(NORMAL);
    if (orig_dest == true)
    {
        contact_normal = PointOfList->GetNormal();
    }
    else
    {
        contact_normal = PointFound->GetNormal();
    }
    
    // Define the discrete contact gap
    Point<3> ProjectedPoint;
    const unsigned int number_nodes = pCond_1->GetGeometry().PointsNumber();
    contact_container.contact_gap.resize(number_nodes);
    contact_container.active_nodes.resize(number_nodes);
    
    for (unsigned int j_index = 0; j_index < number_nodes; j_index++)
    {
        if (orig_dest == true)
        {
            StructuralMechanicsMathUtilities::Project(PointFound->GetPoint(),  pCond_1->GetGeometry()[j_index], ProjectedPoint, contact_container.contact_gap[j_index], contact_normal);
        }
        else
        {
            StructuralMechanicsMathUtilities::Project(PointOfList->GetPoint(), pCond_1->GetGeometry()[j_index], ProjectedPoint, contact_container.contact_gap[j_index], contact_normal);
        }
        
        array_1d<double, 3> result;
        bool inside = pCond_2->GetGeometry().IsInside(ProjectedPoint, result);
        if (inside == true)
        {
            contact_container.active_nodes[j_index] = true;
        }
        else
        {
            contact_container.active_nodes[j_index] = false;
        }
    }
    
    // Define the contact area
    contact_container.contact_area = 0.0;
    
//     KRATOS_WATCH("-----------------------------------------------------------------------------------------------------------------")
//     KRATOS_WATCH(pCond_1->Id());
//     KRATOS_WATCH(pCond_1->GetGeometry());
//     KRATOS_WATCH(pCond_2->Id());
//     KRATOS_WATCH(pCond_2->GetGeometry());
    
    double aux_int = 0.0;
    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = pCond_1->GetGeometry().IntegrationPoints( IntegrationOrder );
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        Point<3> GaussPoint;
        Point<3> GaussPointLocalCoordinates;
        Point<3> ProjectedGaussPoint;
        GaussPointLocalCoordinates.Coordinate(1) = integration_points[PointNumber].X();
        GaussPointLocalCoordinates.Coordinate(2) = integration_points[PointNumber].Y();
        GaussPointLocalCoordinates.Coordinate(3) = integration_points[PointNumber].Z(); // This is supposed to be 0 always, in 1D and 2D
        
        array_1d<double, 3> result;
        GaussPoint = pCond_1->GetGeometry().GlobalCoordinates(result, GaussPointLocalCoordinates);
        
        double dist_aux;
        if (orig_dest == true)
        {
            StructuralMechanicsMathUtilities::Project(PointFound->GetPoint(), GaussPoint,  ProjectedGaussPoint, dist_aux, contact_normal);
        }
        else
        {
            StructuralMechanicsMathUtilities::Project(PointOfList->GetPoint(), GaussPoint, ProjectedGaussPoint, dist_aux, contact_normal);
        }
        
        bool inside = pCond_2->GetGeometry().IsInside(ProjectedGaussPoint, result);
        
//         KRATOS_WATCH(inside);
//         KRATOS_WATCH(result);
//         KRATOS_WATCH(GaussPoint);
//         KRATOS_WATCH(ProjectedGaussPoint);
        
        // Integration weigth
        double IntegrationWeight = integration_points[PointNumber].Weight();
        aux_int += IntegrationWeight;
        
        if (inside == true)
        {
            contact_container.contact_area += IntegrationWeight;
        }
    }
    
    contact_container.contact_area /= aux_int;
    
//     contact_container.print();
        
    ConditionPointers->push_back(contact_container);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ComputeNodesMeanNormal(const bool bidirectional)
{
    // Tolerance
    const double tol = 1.0e-14;
    
    // Initialize normal vectors
    const array_1d<double,3> ZeroNormal = ZeroVector(3);
    
    NodesArrayType& pNodeDestination  = mrDestinationModelPart.Nodes();
    NodesArrayType::iterator it_node_begin = pNodeDestination.ptr_begin();
    NodesArrayType::iterator it_node_end   = pNodeDestination.ptr_end();
    
    for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
    {
        noalias(node_it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
    }
   
    NodesArrayType& pNodeOrigin  = mrOriginModelPart.Nodes();
    if (bidirectional)
    {
        it_node_begin = pNodeOrigin.ptr_begin();
        it_node_end   = pNodeOrigin.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            noalias(node_it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
        }
    }
    
    // Sum all the nodes normals
    ConditionsArrayType& pCondDestination  = mrDestinationModelPart.Conditions();
    ConditionsArrayType::iterator it_cond_begin = pCondDestination.ptr_begin();
    ConditionsArrayType::iterator it_cond_end   = pCondDestination.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_cond_begin; cond_it!=it_cond_end; cond_it++)
    {
        if (cond_it->Is(ACTIVE))
        {
            const array_1d<double, 3> & rNormal = cond_it->GetValue(NORMAL);
            for (unsigned int i = 0; i < cond_it->GetGeometry().PointsNumber(); i++)
            {
                noalias( cond_it->GetGeometry()[i].FastGetSolutionStepValue(NORMAL) ) += rNormal;
            }
        }
    }
    
    if (bidirectional)
    {
        ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
        it_cond_begin = pCondOrigin.ptr_begin();
        it_cond_end   = pCondOrigin.ptr_end();
        
        for(ConditionsArrayType::iterator cond_it = it_cond_begin; cond_it!=it_cond_end; cond_it++)
        {
            if (cond_it->Is(ACTIVE))
            {
                const array_1d<double, 3> & rNormal = cond_it->GetValue(NORMAL);
                for (unsigned int i = 0; i < cond_it->GetGeometry().PointsNumber(); i++)
                {
                    noalias( cond_it->GetGeometry()[i].FastGetSolutionStepValue(NORMAL) ) += rNormal;
                }
            }
        }
    }
    
    // Normalize normal vectors
    it_node_begin = pNodeDestination.ptr_begin();
    it_node_end   = pNodeDestination.ptr_end();
    
    for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
    {
        const double norm = norm_2(node_it->FastGetSolutionStepValue(NORMAL));
        if (norm > tol)
        {
            node_it->FastGetSolutionStepValue(NORMAL)  /= norm;
        }
        // KRATOS_WATCH(Normal);
    }
    
    if (bidirectional)
    {
        it_node_begin = pNodeOrigin.ptr_begin();
        it_node_end   = pNodeOrigin.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            const double norm = norm_2(node_it->FastGetSolutionStepValue(NORMAL));
            if (norm > tol)
            {
                node_it->FastGetSolutionStepValue(NORMAL)  /= norm;
            }
            // KRATOS_WATCH(Normal);
        }
    }
}

/************************************** ACCESS *************************************/
/***********************************************************************************/
    
}  // namespace Kratos.
