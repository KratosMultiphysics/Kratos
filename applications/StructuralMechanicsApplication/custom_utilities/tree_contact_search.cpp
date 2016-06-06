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
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>(mallocation);  
    }
    
    ConditionsArrayType& pCondOrigin = mrOriginModelPart.Conditions();
    it_begin = pCondOrigin.ptr_begin();
    it_end   = pCondOrigin.ptr_end();
    
    for(ConditionsArrayType::iterator cond_it = it_begin; cond_it!=it_end; cond_it++)
    {
        cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>(mallocation); 
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
            cond_it->Set(MASTER, false);
            cond_it->Set(SLAVE, false);   
            
            std::vector<contact_container> * ConditionPointersDestination = cond_it->GetValue(CONTACT_CONTAINERS);
            
            for (unsigned int i =0; i< ConditionPointersDestination->size();i++)
            {
                delete(&((*ConditionPointersDestination)[i].condition));
                
            } 
            
            ConditionPointersDestination->clear();
            
            delete cond_it->GetValue(CONTACT_CONTAINERS); 
            cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>(mallocation); 
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
            cond_it->Set(MASTER, false);
            cond_it->Set(SLAVE, false);   
            
            std::vector<contact_container> * ConditionPointersOrigin = cond_it->GetValue(CONTACT_CONTAINERS);
            
            for (unsigned int i =0; i< ConditionPointersOrigin->size();i++)
            {
                delete(&((*ConditionPointersOrigin)[i].condition));
            } 
            
            ConditionPointersOrigin->clear();
            
            delete cond_it->GetValue(CONTACT_CONTAINERS); 
            cond_it->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>(mallocation); 
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

        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); // TODO: Add the calculation of the normal of the point
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
        
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); // TODO: Add the calculation of the normal of the point
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
        
        PointItem::Pointer pP = PointItem::Pointer(new PointItem(Coord, *(node_it.base()))); // TODO: Add the calculation of the normal of the point
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
        CenterAndRadius(CondOri, Center, Radius, Normal); 
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
        CenterAndRadius(CondOri, Center, Radius, Normal); 
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
        CenterAndRadius(CondOri, Center, Radius, Normal); 
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
                pCondOrigin->Set(MASTER, true);
                pCondOrigin->Set(SLAVE, false); // NOTE: It is supposed to be already false, just in case
            }

            std::vector<contact_container> * ConditionPointersOrigin = pCondOrigin->GetValue(CONTACT_CONTAINERS);
            
//             KRATOS_WATCH(NumberPointsFound); 
            for(unsigned int i = 0; i < NumberPointsFound; i++)
            {   
                Condition::Pointer pCondDestination = PointsFound[i]->GetCondition();
                std::vector<contact_container> * ConditionPointersDestination = pCondDestination->GetValue(CONTACT_CONTAINERS);
                
                if (bidirectional == true)
                {                                        
                    contact_container contact_container_origin;
                    
                    contact_container_origin.condition = &* pCondDestination;
                    
                    // Define the normal to the contact
                    array_1d<double, 3> & contact_normal_orig = pCondOrigin->GetValue(CONTACT_NORMAL);
                    contact_normal_orig = mPointListOrigin[cond_it]->GetNormal();
                    
                    // Define the contact gap
                    Point<3> ProjectedPointDestination;
                    Project(PointsFound[i]->GetPoint(), mPointListOrigin[cond_it]->GetPoint(), ProjectedPointDestination, contact_container_origin.contact_gap, contact_normal_orig);
                    
                    // Define the contact area
                    contact_container_origin.contact_area = 0.0;
                    double aux_int = 0.0;
                    /* Reading integration points */
                    const GeometryType::IntegrationPointsArrayType& integration_points = pCondOrigin->GetGeometry().IntegrationPoints( IntegrationOrder );
                    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
                    {
                        Point<3> GaussPoint;
                        Point<3> GaussPointLocalCoordinates;
                        Point<3> ProjectedGaussPoint;
                        GaussPointLocalCoordinates.Coordinate(1) = integration_points[PointNumber].X();
                        GaussPointLocalCoordinates.Coordinate(2) = integration_points[PointNumber].Y();
                        GaussPointLocalCoordinates.Coordinate(3) = integration_points[PointNumber].Z(); // This is supposed to be 0 always, in 1D and 2D
                        
                        array_1d<double, 3> result;
                        GaussPoint = pCondDestination->GetGeometry().GlobalCoordinates(result, GaussPointLocalCoordinates);
                        
                        double dist_aux;
                        Project(PointsFound[i]->GetPoint(), GaussPoint, ProjectedGaussPoint, dist_aux, contact_normal_orig);
                        
                        bool inside = pCondDestination->GetGeometry().IsInside(ProjectedGaussPoint, result);
                        
                        // Integration weigth
                        double IntegrationWeight = integration_points[PointNumber].Weight();
                        aux_int += IntegrationWeight;
                        
                        if (inside == true)
                        {
                            contact_container_origin.contact_area += IntegrationWeight;
                        }
                    }
                    
                    contact_container_origin.contact_area /= aux_int;
                        
                    ConditionPointersOrigin->push_back(contact_container_origin);
                }
                
                // Set the corresponding flags
                pCondDestination->Set(ACTIVE, true);
                pCondDestination->Set(MASTER, false); // NOTE: It is supposed to be already false, just in case
                pCondDestination->Set(SLAVE, true);
                
                contact_container contact_container_destination;
                    
                contact_container_destination.condition = &* pCondOrigin;
                
                // Define the normal to the contact
                array_1d<double, 3> & contact_normal_dest = pCondDestination->GetValue(CONTACT_NORMAL);
                contact_normal_dest = PointsFound[i]->GetNormal();
                
                // Define the contact gap
                Point<3> ProjectedPointOrigin;
                Project(mPointListOrigin[cond_it]->GetPoint(), PointsFound[i]->GetPoint(), ProjectedPointOrigin, contact_container_destination.contact_gap, contact_normal_dest);
                
                // Define the contact area
                contact_container_destination.contact_area = 0.0;
                double aux_int = 0.0;
                
//                 KRATOS_WATCH("-----------------------------------------------------------------------------------------------------------------")
//                 KRATOS_WATCH(pCondDestination->Id());
//                 KRATOS_WATCH(pCondDestination->GetGeometry());
//                 KRATOS_WATCH(pCondOrigin->Id());
//                 KRATOS_WATCH(pCondOrigin->GetGeometry());
                
                /* Reading integration points */
                const GeometryType::IntegrationPointsArrayType& integration_points = pCondDestination->GetGeometry().IntegrationPoints( IntegrationOrder );
                for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
                {
                    Point<3> GaussPoint;
                    Point<3> GaussPointLocalCoordinates;
                    Point<3> ProjectedGaussPoint;
                    GaussPointLocalCoordinates.Coordinate(1) = integration_points[PointNumber].X();
                    GaussPointLocalCoordinates.Coordinate(2) = integration_points[PointNumber].Y();
                    GaussPointLocalCoordinates.Coordinate(3) = integration_points[PointNumber].Z(); // This is supposed to be 0 always, in 1D and 2D
                    
                    array_1d<double, 3> result;
                    GaussPoint = pCondDestination->GetGeometry().GlobalCoordinates(result, GaussPointLocalCoordinates);
                
                    double dist_aux;
                    Project(mPointListOrigin[cond_it]->GetPoint(), GaussPoint, ProjectedGaussPoint, dist_aux, contact_normal_dest);
                    
                    bool inside = pCondOrigin->GetGeometry().IsInside(ProjectedGaussPoint, result);
                    
//                     KRATOS_WATCH("----------------------------------")
//                     KRATOS_WATCH(inside);
//                     KRATOS_WATCH(result);
//                     KRATOS_WATCH(GaussPoint);
//                     KRATOS_WATCH(ProjectedGaussPoint);
                    
                    // Integration weigth
                    double IntegrationWeight = integration_points[PointNumber].Weight();
                    aux_int += IntegrationWeight;
                    
                    if (inside == true)
                    {
                        contact_container_destination.contact_area += IntegrationWeight;
                    }
//                     KRATOS_WATCH(contact_container_destination.contact_area/2.00); // Divided by two for the line
                }
                
                contact_container_destination.contact_area /= aux_int;
                
                ConditionPointersDestination->push_back(contact_container_destination);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::Project(
        const Point<3>& PointOrigin,
        const Point<3>& PointDestiny,
        Point<3>& PointProjected,
        double & dist,
        const array_1d<double,3> & Normal
        )
{
     array_1d<double,3> vector_points;
     vector_points[0] = PointDestiny.Coordinate(1) - PointOrigin.Coordinate(1);
     vector_points[1] = PointDestiny.Coordinate(2) - PointOrigin.Coordinate(2);
     vector_points[2] = PointDestiny.Coordinate(3) - PointOrigin.Coordinate(3);
    
     dist = vector_points[0] * Normal[0]
          + vector_points[1] * Normal[1]
          + vector_points[2] * Normal[2];
     
     PointProjected.Coordinate(1) = PointDestiny.Coordinate(1) - Normal[0] * dist;
     PointProjected.Coordinate(2) = PointDestiny.Coordinate(2) - Normal[1] * dist;
     PointProjected.Coordinate(3) = PointDestiny.Coordinate(3) - Normal[2] * dist;
}
/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CenterAndRadius(
        const Condition::Pointer pCond,
        Point<3>& Center,
        double& Radius,
        array_1d<double,3> & Normal
        )
{
    Radius = 0.0;
    Center = pCond->GetGeometry().Center();
    Normal = ZeroVector(3);
    
    // TODO: To calculate the normal I am going to use the Newell's method for quadrilateral, I recommend to find some way to compute in a way that it is possible to have the normal in the nodes and use the Nagata Patch
    if (pCond->GetGeometry().PointsNumber() == 2) // A linear line
    {
        array_1d<double,3> v1,v2;

        // Assuming plane X-Y
        v1[0] = 0.0;
        v1[1] = 0.0;
        v1[2] = 1.0;
        
        v2[0] = pCond->GetGeometry()[1].X() - pCond->GetGeometry()[0].X();
        v2[1] = pCond->GetGeometry()[1].Y() - pCond->GetGeometry()[0].Y();
        v2[2] = pCond->GetGeometry()[1].Z() - pCond->GetGeometry()[0].Z();
        
        MathUtils<double>::CrossProduct(Normal,v1,v2);

        double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
        Normal /= NNorm;
    }
    else if (pCond->GetGeometry().PointsNumber() == 3) // A triangle or quadratic line
    {
        if (mdimension == 2) // Quadratic line
        {
            boost::numeric::ublas::bounded_matrix<double, 3, 3 >  matrix_coeficients     = ZeroMatrix(3, 3);
            boost::numeric::ublas::bounded_matrix<double, 3, 3 >  inv_matrix_coeficients = ZeroMatrix(3, 3);
            boost::numeric::ublas::bounded_matrix<double, 3, 1 >  vector_coeficients     = ZeroMatrix(3, 1);
            for (unsigned int i = 0; i < 3; i++)
            {
                matrix_coeficients(i, 0) = pCond->GetGeometry()[i].X() * pCond->GetGeometry()[i].X();
                matrix_coeficients(i, 1) = pCond->GetGeometry()[i].X();
                matrix_coeficients(i, 2) = 1.0;
                
                vector_coeficients(i, 0) = pCond->GetGeometry()[i].Y();
            }
            
            StructuralMechanicsMathUtilities::InvMat3x3(matrix_coeficients, inv_matrix_coeficients);
            
            boost::numeric::ublas::bounded_matrix<double, 3, 1 >  coeficients;
            noalias(coeficients) = prod(inv_matrix_coeficients, vector_coeficients);
            
            Normal[0] =   2.0 * pCond->GetGeometry()[1].X() * coeficients(0, 0) + coeficients(1, 0);
            Normal[1] = - 1.0;
            Normal[2] =   0.0;
            
            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
            Normal /= NNorm;
        }
        else // Triangle
        {
            array_1d<double,3> v1,v2;

            v1[0] = pCond->GetGeometry()[1].X() - pCond->GetGeometry()[0].X();
            v1[1] = pCond->GetGeometry()[1].Y() - pCond->GetGeometry()[0].Y();
            v1[2] = pCond->GetGeometry()[1].Z() - pCond->GetGeometry()[0].Z();

            v2[0] = pCond->GetGeometry()[2].X() - pCond->GetGeometry()[0].X();
            v2[1] = pCond->GetGeometry()[2].Y() - pCond->GetGeometry()[0].Y();
            v2[2] = pCond->GetGeometry()[2].Z() - pCond->GetGeometry()[0].Z();

            MathUtils<double>::CrossProduct(Normal,v1,v2);

            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
            Normal /= NNorm;
        }
    }
    else if (pCond->GetGeometry().PointsNumber() == 4) // A quadrilateral
    {
        // Newell's method
        for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
        {
            unsigned int index_aux = i + 1; 
            if (i == pCond->GetGeometry().PointsNumber() - 1)
            {
                index_aux = 0;
            }
            Normal[0] += (pCond->GetGeometry()[i].Y() - pCond->GetGeometry()[index_aux].Y()) * 
                         (pCond->GetGeometry()[i].Z() - pCond->GetGeometry()[index_aux].Z());
                         
            Normal[1] += (pCond->GetGeometry()[i].Z() - pCond->GetGeometry()[index_aux].Z()) * 
                         (pCond->GetGeometry()[i].X() - pCond->GetGeometry()[index_aux].X());
                         
            Normal[2] += (pCond->GetGeometry()[i].X() - pCond->GetGeometry()[index_aux].X()) * 
                         (pCond->GetGeometry()[i].Y() - pCond->GetGeometry()[index_aux].Y());
        }
        
        double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
        Normal /= NNorm;
    }
    else // The Newell's method can be used, but nodes must be reordered
    {
        KRATOS_THROW_ERROR( std::logic_error, " There is not any method to calculate the normal for this geometry. Number of nodes: ", pCond->GetGeometry().PointsNumber() );
    }
    
    // TODO: Add calculation of radius to geometry.h in june (after release)
    // TODO: Add calculation of normal to geometry.h in june (after release)
    for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
    {
        double dx = Center.Coordinate(1) - pCond->GetGeometry()[i].X();
        double dy = Center.Coordinate(2) - pCond->GetGeometry()[i].Y();
        double dz = Center.Coordinate(3) - pCond->GetGeometry()[i].Z();

        double tmp = dx * dx + dy * dy + dz * dz;

        if(tmp > Radius)
        {
            Radius = tmp;
        }
    }

    Radius = std::sqrt(Radius);
}

/************************************** ACCESS *************************************/
/***********************************************************************************/
    
}  // namespace Kratos.
