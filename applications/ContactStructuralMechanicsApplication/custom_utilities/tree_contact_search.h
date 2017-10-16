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
#include <omp.h>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

/* Custom utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/search_utilities.h"

/* Custom includes*/

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
            "dual_search_check"                    : false,
            "strict_search_check"                  : true,
            "use_exact_integration"                : true
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mAllocationSize = ThisParameters["allocation_size"].GetInt();
        mSearchFactor = ThisParameters["search_factor"].GetDouble();
        mDualSearchCheck = ThisParameters["dual_search_check"].GetBool();
        mStrictSearchCheck = ThisParameters["strict_search_check"].GetBool();
        mUseExactIntegration = ThisParameters["use_exact_integration"].GetBool();
        mSearchTreeType = ConvertSearchTree(ThisParameters["type_search"].GetString());
        mBucketSize = ThisParameters["bucket_size"].GetInt();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->Set(ACTIVE, false);
        }
        
        // Iterate in the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            it_cond->Set(ACTIVE, false);
        }
    }
    
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
    
    void InitializeMortarConditions()
    {
        // Iterate in the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;

            it_cond->GetValue(MAPPING_PAIRS) = boost::shared_ptr<ConditionMap>(new ConditionMap); 
//             it_cond->GetValue(MAPPING_PAIRS)->reserve(mAllocationSize); 
        }
    }
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearScalarMortarConditions()
    {
        ResetContactOperators();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            if (it_node->Is(ACTIVE) == true)
            {
                it_node->Set( ACTIVE, false );
                it_node->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
            }
        }  
    }
    
    void TotalClearComponentsMortarConditions()
    {
        ResetContactOperators();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            if (it_node->Is(ACTIVE) == true)
            {
                it_node->Set( ACTIVE, false );
                noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = ZeroVector(3);
            }
        }  
    }
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void TotalClearALMFrictionlessMortarConditions()
    {        
        ResetContactOperators();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            if (it_node->Is(ACTIVE) == true)
            {
                it_node->Set( ACTIVE, false );
                it_node->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
            }
        }  
    }
    
    /**
     * This function clears the mortar conditions already created (scalar version)
     */
    
    void PartialClearScalarMortarConditions()
    {
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ACTIVE) == false)
            {
                it_node->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
            }
        } 
    }
    
    /**
     * This function clears the mortar conditions already created (components version)
     */
        
    void PartialClearComponentsMortarConditions()
    {
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ACTIVE) == false)
            {
                noalias(it_node->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = ZeroVector(3);
            }
        } 
    }
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void PartialClearALMFrictionlessMortarConditions()
    {
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ACTIVE) == false)
            {
                it_node->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
            }
        } 
    }
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar()
    {
        // Clearing the vector
        mPointListDestination.clear();
        
        // Iterate in the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        // Creating a buffer for parallel vector fill
        const unsigned int num_threads = omp_get_max_threads();
        std::vector<PointVector> points_buffer(num_threads);

        #pragma omp parallel
        {
            const unsigned int thread_id = omp_get_thread_num();

            #pragma omp for
            for(int i = 0; i < num_conditions; i++) 
            {
                auto it_cond = conditions_array.begin() + i;
                
                if (it_cond->Is(MASTER) == true)
                {
                    PointTypePointer p_point = PointTypePointer(new PointItem((*it_cond.base())));
//                     (mPointListDestination).push_back(p_point);
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

    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar()
    {
        const double& delta_time = mrMainModelPart.GetProcessInfo()[DELTA_TIME];
        
        const int num_points = static_cast<int>(mPointListDestination.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_points; i++) 
        {
            mPointListDestination[i]->UpdatePoint(delta_time);
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
        
        const double& delta_time = mrMainModelPart.GetProcessInfo()[DELTA_TIME];
        
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
            for(int i = 0; i < num_conditions; i++) 
            {
                auto it_cond = conditions_array.begin() + i;
                
                if (it_cond->Is(SLAVE) == true)
                {
                    if (mSearchTreeType == KdtreeInRadius)
                    {                  
                        Point<3> center;
                        if (mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X) == true)
                        {
                            center = ContactUtilities::GetHalfJumpCenter(it_cond->GetGeometry(), delta_time); // NOTE: Center in half delta time
                        }
                        else
                        {
                            center = it_cond->GetGeometry().Center(); // NOTE: Real center
                        }
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
                        boost::shared_ptr<ConditionMap>& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
                        Condition::Pointer p_cond_slave = (*it_cond.base()); // MASTER
                        const array_1d<double, 3>& contact_normal_origin = p_cond_slave->GetValue(NORMAL);
                        const double active_check_length = p_cond_slave->GetGeometry().Length() * p_cond_slave->GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
                        
                        // If not active we check if can be potentially in contact
                        if (mUseExactIntegration == false) // LEGACY WAY
                        {
                            for(unsigned int i_pair = 0; i_pair < number_points_found; i_pair++)
                            {   
                                Condition::Pointer p_cond_origin = points_found[i_pair]->GetCondition();
                                
                                const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_origin);
                                
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
                            const unsigned int number_of_nodes = (it_cond->GetGeometry()).size();
                            const unsigned int dimension = (it_cond->GetGeometry()).WorkingSpaceDimension();

                            if (dimension == 2 && number_of_nodes == 2)
                            {
                                MortarKinematicVariables<2> rVariables;
                                MortarOperator<2> rThisMortarConditionMatrices;
                                ExactMortarIntegrationUtility<2, 2> integration_utility = ExactMortarIntegrationUtility<2, 2>(2);
                                
                                for(unsigned int i_pair = 0; i_pair < number_points_found; i_pair++)
                                {   
                                    bool condition_is_active = false;
                                                                
                                    Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                    const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                        
                                    const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master);
                                    
                                    if (condition_checked_right == OK)
                                    {   
                                        condition_is_active = SearchUtilities::CheckExactIntegration<2, 2, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        // If condition is active we add
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                        }
                                        else
                                        {
                                            conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                        }
                                    }
                                    else if (condition_checked_right == AlreadyInTheMap)
                                    {
                                        condition_is_active = SearchUtilities::CheckExactIntegration<2, 2, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->SetActive(p_cond_master);
                                        }
                                    }
                                }
                            }
                            else if (dimension == 3 && number_of_nodes == 3)
                            {
                                MortarKinematicVariables<3> rVariables;
                                MortarOperator<3> rThisMortarConditionMatrices;
                                ExactMortarIntegrationUtility<3, 3> integration_utility = ExactMortarIntegrationUtility<3, 3>(3);
                                
                                for(unsigned int i_pair = 0; i_pair < number_points_found; i_pair++)
                                {   
                                    bool condition_is_active = false;
                                                                
                                    Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                    const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                        
                                    const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master);
                                    
                                    if (condition_checked_right == OK)
                                    {   
                                        condition_is_active = SearchUtilities::CheckExactIntegration<3, 3, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        // If condition is active we add
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                        }
                                        else
                                        {
                                            conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                        }
                                    }
                                    else if (condition_checked_right == AlreadyInTheMap)
                                    {
                                        condition_is_active = SearchUtilities::CheckExactIntegration<3, 3, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->SetActive(p_cond_master);
                                        }
                                    }
                                }
                            }
                            else if (dimension == 3 && number_of_nodes == 4)
                            {
                                MortarKinematicVariables<4> rVariables;
                                MortarOperator<4> rThisMortarConditionMatrices;
                                ExactMortarIntegrationUtility<3, 4> integration_utility = ExactMortarIntegrationUtility<3, 4>(3);
                                
                                for(unsigned int i_pair = 0; i_pair < number_points_found; i_pair++)
                                {   
                                    bool condition_is_active = false;
                                                                
                                    Condition::Pointer p_cond_master = points_found[i_pair]->GetCondition(); // MASTER
                                    const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                                        
                                    const CheckResult condition_checked_right = CheckCondition(conditions_pointers_destination, p_cond_slave, p_cond_master);
                                    
                                    if (condition_checked_right == OK)
                                    {   
                                        condition_is_active = SearchUtilities::CheckExactIntegration<3, 4, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        // If condition is active we add
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->AddNewActiveCondition(p_cond_master);
                                        }
                                        else
                                        {
                                            conditions_pointers_destination->AddNewInactiveCondition(p_cond_master);
                                        }
                                    }
                                    else if (condition_checked_right == AlreadyInTheMap)
                                    {
                                        condition_is_active = SearchUtilities::CheckExactIntegration<3, 4, true>(rVariables, rThisMortarConditionMatrices, integration_utility, p_cond_slave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
                                        
                                        if (condition_is_active == true)
                                        {
                                            conditions_pointers_destination->SetActive(p_cond_master);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED: dimension = " << dimension << " number_of_nodes = " << number_of_nodes << std::endl;
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
    
    /**
     * This function has as pourpose to clean the existing pairs
     */
    
    void CleanMortarConditions()
    {
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            if ( (it_cond)->Is(ACTIVE) == true )
            {
                boost::shared_ptr<ConditionMap>& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
                
                // Initialize geometries
                const array_1d<double, 3>& contact_normal = it_cond->GetValue(NORMAL);
                const double active_check_length = it_cond->GetGeometry().Length() * it_cond->GetProperties().GetValue(ACTIVE_CHECK_FACTOR);
                
                if (mUseExactIntegration == false) // LEGACY WAY
                {
                    for (auto it_pair = conditions_pointers_destination->begin(); it_pair != conditions_pointers_destination->end(); ++it_pair )
                    {
                        SearchUtilities::ContactContainerFiller<false>(conditions_pointers_destination, (*it_cond.base()), (it_pair->first), contact_normal, (it_pair->first)->GetValue(NORMAL), active_check_length, mDualSearchCheck, mStrictSearchCheck);
                    }
                }
                else
                {
                    const unsigned int number_of_nodes = (it_cond->GetGeometry()).size();
                    const unsigned int dimension = (it_cond->GetGeometry()).WorkingSpaceDimension();
                 
                    const array_1d<double, 3>& contact_normal_origin = it_cond->GetValue(NORMAL);
                    
                    if (dimension == 2 && number_of_nodes == 2)
                    {
                        SearchUtilities::ExactContactContainerChecker<2,2>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                    }
                    else if (dimension == 3 && number_of_nodes == 3)
                    {
                        SearchUtilities::ExactContactContainerChecker<3,3>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                    }
                    else if (dimension == 3 && number_of_nodes == 4)
                    {
                        SearchUtilities::ExactContactContainerChecker<3,4>(conditions_pointers_destination, it_cond->GetGeometry(), contact_normal_origin, active_check_length); 
                    }
                    else
                    {
                        KRATOS_ERROR << "INTEGRATION NOT IMPLEMENTED: dimension = " << dimension << " number_of_nodes = " << number_of_nodes << std::endl;
                    }
                }
            }
        }
    }
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions()
    {
        // Iterate in the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(SLAVE) == true && it_cond->Is(ACTIVE) == true)
            {
                KRATOS_WATCH(it_cond->Id());
                KRATOS_WATCH(it_cond->GetGeometry());
                
                boost::shared_ptr<ConditionMap>& conditions_pointers_destination = it_cond->GetValue(MAPPING_PAIRS);
                KRATOS_WATCH(conditions_pointers_destination->size());
                conditions_pointers_destination->print();
            }
        }
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            if (it_node->Is(ACTIVE) == true)
            {
                std::cout << "Node: " << it_node->Id() << " is active" << std::endl;
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
    
    static inline CheckResult CheckCondition(
        boost::shared_ptr<ConditionMap>& ConditionPointers1,
        const Condition::Pointer& pCond1,
        const Condition::Pointer& pCond2
        )
    {
        if (pCond1 == pCond2) // Avoiding "auto self-contact"
        {
            return Fail;
        }
        
        
        if (((pCond1->GetValue(ELEMENT_POINTER) != nullptr) && (pCond2->GetValue(ELEMENT_POINTER) != nullptr)) == true)
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

        if (pCond2->Is(SLAVE) == true) // Otherwise will not be necessary to check
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
    
    /**  
     * Calculates the minimal distance between one node and its center 
     * @return The radius of the geometry 
     */ 
    
    static inline double Radius(GeometryType& ThisGeometry) 
    { 
        double radius = 0.0; 
        const Point<3>& center = ThisGeometry.Center(); 
         
        for(unsigned int i_node = 0; i_node < ThisGeometry.PointsNumber(); i_node++) 
        { 
            const array_1d<double, 3> aux_vector = center.Coordinates() - ThisGeometry[i_node].Coordinates();
             
            const double aux_value = inner_prod(aux_vector, aux_vector); 
 
            if(aux_value > radius) 
            { 
                radius = aux_value; 
            } 
        } 
 
        return std::sqrt(radius); 
    } 

    /**
     * This resets the contact operators
     */
        
    void ResetContactOperators()
    {
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            if (it_cond->Is(SLAVE) == true && it_cond->Is(ACTIVE) == true)
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
    bool mDualSearchCheck;                    // The search is done reciprocally
    bool mStrictSearchCheck;                  // The search is done requiring IsInside as true
    bool mUseExactIntegration;                // The search filter the results with the exact integration
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
