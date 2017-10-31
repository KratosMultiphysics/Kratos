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

// System includes

// External includes

// Project includes

/* Custom utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/tree_contact_search.h"

namespace Kratos
{
TreeContactSearch::TreeContactSearch(
        ModelPart & rMainModelPart,
        Parameters ThisParameters
        ):mrMainModelPart(rMainModelPart.GetSubModelPart("Contact")),
        mDimension(rMainModelPart.GetProcessInfo()[DOMAIN_SIZE])
{        
    Parameters DefaultParameters = Parameters(R"(
    {
        "allocation_size"                      : 1000, 
        "bucket_size"                          : 4, 
        "search_factor"                        : 2.0, 
        "type_search"                          : "InRadius", 
        "dual_search_check"                    : false,
        "strict_search_check"                  : true,
        "use_exact_integration"                : true,
        "inverted_search"                      : false
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

    mAllocationSize = ThisParameters["allocation_size"].GetInt();
    mSearchFactor = ThisParameters["search_factor"].GetDouble();
    mDualSearchCheck = ThisParameters["dual_search_check"].GetBool();
    mStrictSearchCheck = ThisParameters["strict_search_check"].GetBool();
    mUseExactIntegration = ThisParameters["use_exact_integration"].GetBool();
    mInvertedSearch = ThisParameters["inverted_search"].GetBool();
    mSearchTreeType = ConvertSearchTree(ThisParameters["type_search"].GetString());
    mBucketSize = ThisParameters["bucket_size"].GetInt();
    
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->Set(ACTIVE, false);
    }
    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        it_cond->Set(ACTIVE, false);
    }
}   
    
/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InitializeMortarConditions()
{
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;

        if (it_cond->Has(MAPPING_PAIRS) == false) it_cond->SetValue(MAPPING_PAIRS, ConditionMap::Pointer(new ConditionMap)); 
//             it_cond->GetValue(MAPPING_PAIRS)->reserve(mAllocationSize); 
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::TotalClearScalarMortarConditions()
{
    TotalResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
    }  
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::TotalClearComponentsMortarConditions()
{
    TotalResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = ZeroVector(3);
    }  
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::TotalClearALMFrictionlessMortarConditions()
{        
    TotalResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
    }  
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::PartialClearScalarMortarConditions()
{
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
    } 
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::PartialClearComponentsMortarConditions()
{
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = ZeroVector(3);
    } 
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::PartialClearALMFrictionlessMortarConditions()
{
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
    } 
}
    
/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CreatePointListMortar()
{
    // Clearing the vector
    mPointListDestination.clear();
    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<PointVector> points_buffer(num_threads);

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        #pragma omp for
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(MASTER) == !mInvertedSearch)
            {
                const PointTypePointer& p_point = PointTypePointer(new PointItem((*it_cond.base())));
                (points_buffer[thread_id]).push_back(p_point);
            }
        }
        
        // Combine buffers together
        #pragma omp single
        {
            for( auto& point_buffer : points_buffer)
            {
                std::move(point_buffer.begin(),point_buffer.end(),back_inserter(mPointListDestination));
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdatePointListMortar()
{
    const double& delta_time = mrMainModelPart.GetProcessInfo()[DELTA_TIME];
    
    const int num_points = static_cast<int>(mPointListDestination.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_points; ++i) mPointListDestination[i]->UpdatePoint(delta_time);
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::UpdateMortarConditions()
{        
    // We update the list of points
    UpdatePointListMortar();
    
    // Calculate the mean of the normal in all the nodes
    ContactUtilities::ComputeNodesMeanNormalModelPart(mrMainModelPart); 
    
    const double& delta_time = mrMainModelPart.GetProcessInfo()[DELTA_TIME];
    
    // We check if we are in a dynamic or static case
    const bool dynamic = mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X) ;
    
    // Taking the ACTIVE_CHECK_FACTOR
    const double& active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    
//         #pragma omp parallel 
//         {
        // Initialize values
        PointVector points_found(mAllocationSize);
        unsigned int number_points_found = 0;    
        
        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTree tree_points(mPointListDestination.begin(), mPointListDestination.end(), mBucketSize);
        
        // Iterate in the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

//             #pragma omp for 
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(SLAVE) == !mInvertedSearch)
            {
                if (mSearchTreeType == KdtreeInRadius)
                {
                    GeometryType& geometry = it_cond->GetGeometry();
                    const Point& center = dynamic ? ContactUtilities::GetHalfJumpCenter(geometry, delta_time) : geometry.Center(); // NOTE: Center in half delta time or real center
                    
                    const double search_radius = mSearchFactor * Radius(it_cond->GetGeometry());

                    number_points_found = tree_points.SearchInRadius(center, search_radius, points_found.begin(), mAllocationSize);
                }
                else if (mSearchTreeType == KdtreeInBox)
                {
                    // Auxiliar values
                    const double length_search = mSearchFactor * it_cond->GetGeometry().Length();
                    
                    // Compute max/min points
                    Node<3> min_point, max_point;
                    it_cond->GetGeometry().BoundingBox(min_point, max_point);
                    
                    // Get the normal in the extrema points
                    Vector N_min, N_max;
                    GeometryType::CoordinatesArrayType local_point_min, local_point_max;
                    it_cond->GetGeometry().PointLocalCoordinates( local_point_min, min_point.Coordinates( ) ) ;
                    it_cond->GetGeometry().PointLocalCoordinates( local_point_max, max_point.Coordinates( ) ) ;
                    it_cond->GetGeometry().ShapeFunctionsValues( N_min, local_point_min );
                    it_cond->GetGeometry().ShapeFunctionsValues( N_max, local_point_max );
                
                    const array_1d<double,3> normal_min = MortarUtilities::GaussPointUnitNormal(N_min, it_cond->GetGeometry());
                    const array_1d<double,3> normal_max = MortarUtilities::GaussPointUnitNormal(N_max, it_cond->GetGeometry());
                    
                    ContactUtilities::ScaleNode<Node<3>>(min_point, normal_min, length_search);
                    ContactUtilities::ScaleNode<Node<3>>(max_point, normal_max, length_search);
                    
                    number_points_found = tree_points.SearchInBox(min_point, max_point, points_found.begin(), mAllocationSize);
                }
                else
                {
                    KRATOS_ERROR << " The type search is not implemented yet does not exist!!!!. SearchTreeType = " << mSearchTreeType << std::endl;
                }
                
                if (number_points_found > 0)
                {                           
                    ConditionMap::Pointer& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
                    Condition::Pointer p_cond_slave = (*it_cond.base()); // MASTER
                    const array_1d<double, 3>& contact_normal_origin = p_cond_slave->GetValue(NORMAL);
                    const GeometryType& this_geometry = p_cond_slave->GetGeometry();
                    const double active_check_length = this_geometry.Length() *active_check_factor;
                    
                    // If not active we check if can be potentially in contact
                    if (mUseExactIntegration == false) // LEGACY WAY
                    {
                        for(unsigned int i_pair = 0; i_pair < number_points_found; ++i_pair)
                        {   
                            Condition::Pointer p_cond_origin = points_found[i_pair]->GetCondition();
                            
                            const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_origin, mInvertedSearch);
                            
                            if (condition_checked_right == OK)
                            {    
                                SearchUtilities::ContactContainerFiller<true>(conditions_pointers_destination, p_cond_slave, p_cond_origin, contact_normal_origin, p_cond_origin->GetValue(NORMAL), active_check_length, mDualSearchCheck, mStrictSearchCheck); 
                            }
                            else if (condition_checked_right == AlreadyInTheMap)
                            {
                                SearchUtilities::ContactContainerFiller<false>(conditions_pointers_destination, p_cond_slave, p_cond_origin, contact_normal_origin, p_cond_origin->GetValue(NORMAL), active_check_length, mDualSearchCheck, mStrictSearchCheck); 
                            }
                        }
                    }
                    else
                    {
                        if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
                        {
                            MortarKinematicVariables<2> rVariables;
                            MortarOperator<2> rThisMortarConditionMatrices;
                            ExactMortarIntegrationUtility<2, 2> integration_utility = ExactMortarIntegrationUtility<2, 2>(2);
                            
                            for(unsigned int i_pair = 0; i_pair < number_points_found; ++i_pair)
                            {   
                                bool condition_is_active = false;
                                                            
                                Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                    
                                const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master, mInvertedSearch);
                                
                                if (condition_checked_right == OK)
                                {   
                                    condition_is_active = SearchUtilities::CheckExactIntegration<2, 2, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    // If condition is active we add
                                    if (condition_is_active) conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                    else conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                }
                                else if (condition_checked_right == AlreadyInTheMap)
                                {
                                    condition_is_active = SearchUtilities::CheckExactIntegration<2, 2, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    if (condition_is_active) conditions_pointers_destination->SetActive(p_cond_master);
                                }
                            }
                        }
                        else if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
                        {
                            MortarKinematicVariables<3> rVariables;
                            MortarOperator<3> rThisMortarConditionMatrices;
                            ExactMortarIntegrationUtility<3, 3> integration_utility = ExactMortarIntegrationUtility<3, 3>(3);
                            
                            for(unsigned int i_pair = 0; i_pair < number_points_found; ++i_pair)
                            {   
                                bool condition_is_active = false;
                                                            
                                Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                    
                                const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master, mInvertedSearch);
                                
                                if (condition_checked_right == OK)
                                {   
                                    condition_is_active = SearchUtilities::CheckExactIntegration<3, 3, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    // If condition is active we add
                                    if (condition_is_active == true) conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                    else conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                }
                                else if (condition_checked_right == AlreadyInTheMap)
                                {
                                    condition_is_active = SearchUtilities::CheckExactIntegration<3, 3, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    if (condition_is_active == true) conditions_pointers_destination->SetActive(p_cond_master);
                                }
                            }
                        }
                        else if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
                        {
                            MortarKinematicVariables<4> rVariables;
                            MortarOperator<4> rThisMortarConditionMatrices;
                            ExactMortarIntegrationUtility<3, 4> integration_utility = ExactMortarIntegrationUtility<3, 4>(3);
                            
                            for(unsigned int i_pair = 0; i_pair < number_points_found; ++i_pair)
                            {   
                                bool condition_is_active = false;
                                                            
                                Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                    
                                const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master, mInvertedSearch);
                                
                                if (condition_checked_right == OK)
                                {   
                                    condition_is_active = SearchUtilities::CheckExactIntegration<3, 4, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    // If condition is active we add
                                    if (condition_is_active) conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                    else conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                }
                                else if (condition_checked_right == AlreadyInTheMap)
                                {
                                    condition_is_active = SearchUtilities::CheckExactIntegration<3, 4, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                    
                                    if (condition_is_active) conditions_pointers_destination->SetActive(p_cond_master);
                                }
                            }
                        }
                        else
                        {
                            KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED: dimension = " << this_geometry.WorkingSpaceDimension() << " number_of_nodes = " << this_geometry.size() << std::endl;
                        }
                    }
                
                    if ((conditions_pointers_destination->size() > 0) && 
                        (conditions_pointers_destination->AtLeastOnePairActive() == true))
                    {                        
                        it_cond->Set(ACTIVE, true);
                    }
                }
            }
        }
//         }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CleanMortarConditions()
{
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    const double& active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    
    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        const GeometryType& this_geometry = it_cond->GetGeometry();
        if ( (it_cond)->Is(ACTIVE) == true )
        {
            ConditionMap::Pointer& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
            
            // Initialize geometries
            const array_1d<double, 3>& contact_normal = it_cond->GetValue(NORMAL);
            const double active_check_length = this_geometry.Length() * active_check_factor;
            
            if (mUseExactIntegration == false) // LEGACY WAY
            {
                for (auto it_pair = conditions_pointers_destination->begin(); it_pair != conditions_pointers_destination->end(); ++it_pair )
                {
                    SearchUtilities::ContactContainerFiller<false>(conditions_pointers_destination, (*it_cond.base()), (it_pair->first), contact_normal, (it_pair->first)->GetValue(NORMAL), active_check_length, mDualSearchCheck, mStrictSearchCheck);
                }
            }
            else
            {
                const array_1d<double, 3>& contact_normal_origin = it_cond->GetValue(NORMAL);
                
                if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
                {
                    SearchUtilities::ExactContactContainerChecker<2,2>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                }
                else if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
                {
                    SearchUtilities::ExactContactContainerChecker<3,3>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                }
                else if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
                {
                    SearchUtilities::ExactContactContainerChecker<3,4>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                }
                else
                {
                    KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED: dimension = " << this_geometry.WorkingSpaceDimension() << " number_of_nodes = " << this_geometry.size() << std::endl;
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::CheckMortarConditions()
{
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        
        ConditionMap::Pointer& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
        if (conditions_pointers_destination->size() > 0) KRATOS_WATCH(conditions_pointers_destination->size());
        
        if (it_cond->Is(SLAVE) == true && it_cond->Is(ACTIVE) == true)
        {
            KRATOS_WATCH(it_cond->Id());
            KRATOS_WATCH(it_cond->GetGeometry());
            conditions_pointers_destination->print();
        }
    }
    
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());
    
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
        if (it_node->Is(ACTIVE) == true) std::cout << "Node: " << it_node->Id() << " is active" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::InvertSearch()
{
    mInvertedSearch = !mInvertedSearch;
}

/***********************************************************************************/
/***********************************************************************************/

inline CheckResult TreeContactSearch::CheckCondition(
    ConditionMap::Pointer& ConditionPointers1,
    const Condition::Pointer& pCond1,
    const Condition::Pointer& pCond2,
    const bool InvertedSearch
    )
{
    if (pCond1 == pCond2) // Avoiding "auto self-contact"
    {
        return Fail;
    }
    
    if (((pCond1->Has(ELEMENT_POINTER)) && (pCond2->Has(ELEMENT_POINTER))) == true)
    {
        if ((pCond1->GetValue(ELEMENT_POINTER) != pCond2->GetValue(ELEMENT_POINTER)) == false) // Avoiding "auto element contact"
        {
            return Fail;
        }
    }

    // Avoid conditions oriented in the same direction
    const double tolerance = 1.0e-16;
    if (norm_2(pCond1->GetValue(NORMAL) - pCond2->GetValue(NORMAL)) < tolerance)
    {
        return Fail;
    }

    if (pCond2->Is(SLAVE) == !InvertedSearch) // Otherwise will not be necessary to check
    {
        auto& condition_pointers2 = pCond2->GetValue(MAPPING_PAIRS);
        
        if (condition_pointers2->find(pCond1) != condition_pointers2->end())
        {
            return Fail;
        }
    }
    
    // To avoid to repeat twice the same condition 
    if (ConditionPointers1->find(pCond2) != ConditionPointers1->end())
    {
        return AlreadyInTheMap;
    }

    return OK;
}

/***********************************************************************************/
/***********************************************************************************/

inline double TreeContactSearch::Radius(GeometryType& ThisGeometry) 
{ 
    double radius = 0.0; 
    const Point& center = ThisGeometry.Center(); 
        
    for(unsigned int i_node = 0; i_node < ThisGeometry.PointsNumber(); ++i_node) 
    { 
        const array_1d<double, 3>& aux_vector = center.Coordinates() - ThisGeometry[i_node].Coordinates();
            
        const double aux_value = inner_prod(aux_vector, aux_vector); 

        if(aux_value > radius) radius = aux_value; 
    } 

    return std::sqrt(radius); 
} 

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::ResetContactOperators()
{
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        if (it_cond->Is(SLAVE) == !mInvertedSearch && it_cond->Is(ACTIVE) == true)
        {
            it_cond->Set(ACTIVE, false);
            
            auto& condition_pointers = it_cond->GetValue(MAPPING_PAIRS);
            
            if (condition_pointers != nullptr)
            {
                condition_pointers->clear();
//                     condition_pointers->reserve(mAllocationSize); 
            }
        }
    }   
}

/***********************************************************************************/
/***********************************************************************************/

void TreeContactSearch::TotalResetContactOperators()
{
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());
    
    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        it_cond->Set(ACTIVE, false);
        auto& condition_pointers = it_cond->GetValue(MAPPING_PAIRS);
        if (condition_pointers != nullptr)  condition_pointers->clear();
    }   
}

/***********************************************************************************/
/***********************************************************************************/

SearchTreeType TreeContactSearch::ConvertSearchTree(const std::string& str)
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
}  // namespace Kratos.
