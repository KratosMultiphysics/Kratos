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
#include "utilities/mortar_utilities.h"

/* Custom utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/tree_contact_search.h"

namespace Kratos
{
template<unsigned int TDim, unsigned int TNumNodes>
TreeContactSearch<TDim, TNumNodes>::TreeContactSearch( 
    ModelPart & rMainModelPart, 
    Parameters ThisParameters 
    ):mrMainModelPart(rMainModelPart), 
      mThisParameters(ThisParameters)
{        
    KRATOS_ERROR_IF(mrMainModelPart.HasSubModelPart("Contact") == false) << "WARNING:: Please add the Contact submodelpart to your modelpart list" << std::endl;
    
    Parameters DefaultParameters = Parameters(R"(
    {
        "allocation_size"                      : 1000, 
        "bucket_size"                          : 4, 
        "search_factor"                        : 2.0, 
        "type_search"                          : "InRadius", 
        "check_gap"                            : "MappingCheck", 
        "condition_name"                       : "",  
        "final_string"                         : "",  
        "inverted_search"                      : false
    })" );
    
    mThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mCheckGap = ConvertCheckGap(mThisParameters["check_gap"].GetString());
    mInvertedSearch = mThisParameters["inverted_search"].GetBool();
    
    // Check if the computing contact submodelpart
    if (mrMainModelPart.HasSubModelPart("ComputingContact") == false) // We check if the submodelpart where the actual conditions used to compute contact are going to be computed 
        mrMainModelPart.CreateSubModelPart("ComputingContact");
    else { // We clean the existing modelpart
        ModelPart& computing_rcontact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact"); 
        ConditionsArrayType& conditions_array = computing_rcontact_model_part.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; ++i) 
            (conditions_array.begin() + i)->Set(TO_ERASE, true);
        
        mrMainModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    }
    
    // Updating the base condition
    mConditionName = mThisParameters["condition_name"].GetString();
    if (mConditionName == "") 
        mCreateAuxiliarConditions = false;
    else {
        mCreateAuxiliarConditions = true;
        mConditionName.append("Condition"); 
        mConditionName.append(std::to_string(TDim));
        mConditionName.append("D"); 
        mConditionName.append(std::to_string(TNumNodes)); 
        mConditionName.append("N"); 
        mConditionName.append(mThisParameters["final_string"].GetString());
    }
    
    // We iterate over the nodes
    NodesArrayType& nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->Set(ACTIVE, false);
        if (mCheckGap == MappingCheck) it_node->SetValue(NORMAL_GAP, 0.0);
    }
    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.GetSubModelPart("Contact").Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) 
        (conditions_array.begin() + i)->Set(ACTIVE, false);
}   
    
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::InitializeMortarConditions()
{    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.GetSubModelPart("Contact").Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

    #pragma omp parallel for 
    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = conditions_array.begin() + i;

        if (it_cond->Has(INDEX_SET) == false) it_cond->SetValue(INDEX_SET, boost::make_shared<IndexSet>()); 
//             it_cond->GetValue(INDEX_SET)->reserve(mThisParameters["allocation_size"].GetInt()); 
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::ClearScalarMortarConditions()
{
    ResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
        (nodes_array.begin() + i)->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::ClearComponentsMortarConditions()
{
    ResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
        noalias((nodes_array.begin() + i)->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = ZeroVector(3);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::ClearALMFrictionlessMortarConditions()
{        
    ResetContactOperators();
    
    NodesArrayType& nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
        (nodes_array.begin() + i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0;
}
    
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::CreatePointListMortar()
{
    // Clearing the vector
    mPointListDestination.clear();
    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.GetSubModelPart("Contact").Conditions();

    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<PointVector> points_buffer(num_threads);

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        #pragma omp for
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(MASTER) == !mInvertedSearch) {
                const PointTypePointer& p_point = PointTypePointer(new PointItem((*it_cond.base())));
                (points_buffer[thread_id]).push_back(p_point);
            }
        }
        
        // Combine buffers together
        #pragma omp single
        {
            for( auto& point_buffer : points_buffer)
                std::move(point_buffer.begin(),point_buffer.end(),back_inserter(mPointListDestination));
        }
    }
    
#ifdef KRATOS_DEBUG
    // NOTE: We check the list
    for (unsigned int i_point = 0; i_point < mPointListDestination.size(); ++i_point )
        mPointListDestination[i_point]->Check();
//     std::cout << "The list is properly built" << std::endl;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::UpdatePointListMortar()
{
    // We initialize the acceleration
    if (mrMainModelPart.GetProcessInfo()[STEP] == 1) 
        InitializeAcceleration();
    
    // We check if we are in a dynamic or static case
    const bool dynamic = mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X); 
    const double delta_time = (dynamic) ? mrMainModelPart.GetProcessInfo()[DELTA_TIME] : 0.0;
    
    // We compute the delta displacement
    ContactUtilities::ComputeStepJump(mrMainModelPart.GetSubModelPart("Contact"), delta_time);
    
    if (mCheckGap == MappingCheck && dynamic == true) {
        NodesArrayType& update_nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(update_nodes_array.size()); ++i) 
            noalias((update_nodes_array.begin() + i)->Coordinates()) += (update_nodes_array.begin() + i)->GetValue(DELTA_COORDINATES);
    }
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(mPointListDestination.size()); ++i) mPointListDestination[i]->UpdatePoint();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::UpdateMortarConditions()
{            
    // We update the list of points
    UpdatePointListMortar();
    
    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mrMainModelPart.GetSubModelPart("Contact")); 
    
    // We get the computing model part
    std::size_t condition_id = ReorderConditionsIds();
    ModelPart& computing_rcontact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact"); 
    
    // We check if we are in a dynamic or static case
    const bool dynamic = mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X);
    
    // Some auxiliar values
    const unsigned int allocation_size = mThisParameters["allocation_size"].GetInt();           // Allocation size for the vectors and max number of potential results 
    const double search_factor = mThisParameters["search_factor"].GetDouble();                  // The search factor to be considered 
    SearchTreeType type_search = ConvertSearchTree(mThisParameters["type_search"].GetString()); // The search tree considered
    unsigned int bucket_size = mThisParameters["bucket_size"].GetInt();                         // Bucket size for kd-tree
    
    // Create a tree
    // It will use a copy of mNodeList (a std::vector which contains pointers)
    // Copying the list is required because the tree will reorder it for efficiency
    KDTree tree_points(mPointListDestination.begin(), mPointListDestination.end(), bucket_size);
    
    // Iterate in the conditions
    ModelPart& rcontact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ConditionsArrayType& conditions_array = rcontact_model_part.Conditions();
    const int num_conditions = static_cast<int>(conditions_array.size());

//     #pragma omp parallel for firstprivate(tree_points) // FIXME: Make me parallel!!! 
    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = conditions_array.begin() + i;
        
        if (it_cond->Is(SLAVE) == !mInvertedSearch) {
            // Initialize values
            PointVector points_found(allocation_size);
            
            unsigned int number_points_found = 0;    
            
            if (type_search == KdtreeInRadius) {
                GeometryType& geometry = it_cond->GetGeometry();
                const Point& center = dynamic ? ContactUtilities::GetHalfJumpCenter(geometry) : geometry.Center(); // NOTE: Center in half delta time or real center
                
                const double search_radius = search_factor * Radius(it_cond->GetGeometry());

                number_points_found = tree_points.SearchInRadius(center, search_radius, points_found.begin(), allocation_size);
            }
            else if (type_search == KdtreeInBox) {
                // Auxiliar values
                const double length_search = search_factor * it_cond->GetGeometry().Length();
                
                // Compute max/min points
                NodeType min_point, max_point;
                it_cond->GetGeometry().BoundingBox(min_point, max_point);
                
                // Get the normal in the extrema points
                Vector N_min, N_max;
                GeometryType::CoordinatesArrayType local_point_min, local_point_max;
                it_cond->GetGeometry().PointLocalCoordinates( local_point_min, min_point.Coordinates( ) ) ;
                it_cond->GetGeometry().PointLocalCoordinates( local_point_max, max_point.Coordinates( ) ) ;
                it_cond->GetGeometry().ShapeFunctionsValues( N_min, local_point_min );
                it_cond->GetGeometry().ShapeFunctionsValues( N_max, local_point_max );
            
                const array_1d<double,3>& normal_min = MortarUtilities::GaussPointUnitNormal(N_min, it_cond->GetGeometry());
                const array_1d<double,3>& normal_max = MortarUtilities::GaussPointUnitNormal(N_max, it_cond->GetGeometry());
                
                ContactUtilities::ScaleNode<NodeType>(min_point, normal_min, length_search);
                ContactUtilities::ScaleNode<NodeType>(max_point, normal_max, length_search);
                
                number_points_found = tree_points.SearchInBox(min_point, max_point, points_found.begin(), allocation_size);
            }
            else
                KRATOS_ERROR << " The type search is not implemented yet does not exist!!!!. SearchTreeType = " << type_search << std::endl;
            
            if (number_points_found > 0) {   
                // We resize the vector to the actual size
//                 points_found.resize(number_points_found); // NOTE: May be ineficient
                
            #ifdef KRATOS_DEBUG
                // NOTE: We check the list
                for (unsigned int i_point = 0; i_point < number_points_found; ++i_point )
                    points_found[i_point]->Check();
//                 std::cout << "The search is properly done" << std::endl;
            #endif
                
                IndexSet::Pointer indexes_set = it_cond->GetValue(INDEX_SET);
                
                // If not active we check if can be potentially in contact
                if (mCheckGap == MappingCheck) {
                    for (unsigned int i_point = 0; i_point < number_points_found; ++i_point ) {
                        Condition::Pointer p_cond_master = points_found[i_point]->GetCondition();
                        const CheckResult condition_checked_right = CheckCondition(indexes_set, (*it_cond.base()), p_cond_master, mInvertedSearch);
                        if (condition_checked_right == OK) indexes_set->AddId(p_cond_master->Id());
                    }
                }
                else
                    AddPotentialPairing(computing_rcontact_model_part, condition_id, (*it_cond.base()), points_found, number_points_found, indexes_set);
            }
        }
    }

    // We map the Coordinates to the slave side from the master
    if (mCheckGap == MappingCheck) CheckPairing(computing_rcontact_model_part, condition_id);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::AddPairing(
    ModelPart& rComputingModelPart,
    std::size_t& ConditionId,
    Condition::Pointer pCondSlave,
    Condition::Pointer pCondMaster
    )
{    
    if (mCreateAuxiliarConditions == true) { // We add the ID and we create a new auxiliar condition
        ++ConditionId;
        Condition::Pointer p_auxiliar_condition = rComputingModelPart.CreateNewCondition(mConditionName, ConditionId, pCondSlave->GetGeometry(), pCondSlave->pGetProperties());
        // We set the geometrical values
        p_auxiliar_condition->SetValue(PAIRED_GEOMETRY, pCondMaster->pGetGeometry());
        p_auxiliar_condition->SetValue(NORMAL, pCondSlave->GetValue(NORMAL));
        p_auxiliar_condition->SetValue(PAIRED_NORMAL, pCondMaster->GetValue(NORMAL));
        // We activate the condition and initialize it
        p_auxiliar_condition->Set(ACTIVE, true);
        p_auxiliar_condition->Initialize();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::AddPairing(
    ModelPart& rComputingModelPart,
    std::size_t& ConditionId,
    Condition::Pointer pCondSlave,
    Condition::Pointer pCondMaster,
    IndexSet::Pointer IndexesSet
    )
{
    IndexesSet->AddId(pCondMaster->Id());
    
    AddPairing(rComputingModelPart, ConditionId, pCondSlave, pCondMaster);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::CheckMortarConditions()
{
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.GetSubModelPart("Contact").Conditions();

    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
        auto it_cond = conditions_array.begin() + i;
        
        if (it_cond->Has(INDEX_SET) == true) {
            IndexSet::Pointer ids_destination = it_cond->GetValue(INDEX_SET);
            if (ids_destination->size() > 0) {
                std::cout << "Origin condition ID:" << it_cond->Id() << " Number of pairs: " << ids_destination->size() << std::endl;
                std::cout << ids_destination->Info();
            }
        }
    }
    
    NodesArrayType& nodes_array = mrMainModelPart.GetSubModelPart("Contact").Nodes();
    
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        if (it_node->Is(ACTIVE) == true) std::cout << "Node: " << it_node->Id() << " is active" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::InvertSearch()
{
    mInvertedSearch = !mInvertedSearch;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline double TreeContactSearch<TDim, TNumNodes>::GetMaxNodalH()
{
    // We iterate over the nodes
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    
    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<double> max_vector(num_threads, 0.0);
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        const double nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);

        const int id = OpenMPUtils::ThisThread();
        
        if (nodal_h > max_vector[id]) 
            max_vector[id] = nodal_h;
    }
    
    return *std::max_element(max_vector.begin(), max_vector.end());
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline double TreeContactSearch<TDim, TNumNodes>::GetMeanNodalH()
{
    // We iterate over the nodes
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();
    
    double sum_nodal_h = 0.0;
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        const double nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);

        #pragma omp atomic
        sum_nodal_h += nodal_h;
    }
    
    return sum_nodal_h/static_cast<double>(nodes_array.size());
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline void TreeContactSearch<TDim, TNumNodes>::InitializeAcceleration()
{
    // We iterate over the nodes
    NodesArrayType& nodes_array = mrMainModelPart.Nodes();

    // We check if we are in a dynamic or static case
    const bool dynamic = mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X);
    if (dynamic == true) {
        array_1d<double, 3> volume_acceleration = ZeroVector(3);
        auto& r_props = mrMainModelPart.rProperties();
        for (auto& prop : r_props) {
            if (prop.Has(VOLUME_ACCELERATION) == true) {
                volume_acceleration = prop.GetValue(VOLUME_ACCELERATION);
                if (norm_2(volume_acceleration) > 0.0) break;
            }
        }
        
        // Initialize the acceleration when VOLUME_ACCELERATION to be able to predict position in the first step
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            const array_1d<double, 3>& nodal_volume_acceleration = it_node->FastGetSolutionStepValue(VOLUME_ACCELERATION);
            if (norm_2(volume_acceleration) > 0.0) {
                if (!it_node->IsFixed(ACCELERATION_X)) it_node->FastGetSolutionStepValue(ACCELERATION_X) = volume_acceleration[0];
                if (!it_node->IsFixed(ACCELERATION_Y)) it_node->FastGetSolutionStepValue(ACCELERATION_Y) = volume_acceleration[1];
                if (!it_node->IsFixed(ACCELERATION_Z)) it_node->FastGetSolutionStepValue(ACCELERATION_Z) = volume_acceleration[2];
            }
            else if (norm_2(nodal_volume_acceleration) > 0.0) {
                if (!it_node->IsFixed(ACCELERATION_X)) it_node->FastGetSolutionStepValue(ACCELERATION_X) = nodal_volume_acceleration[0];
                if (!it_node->IsFixed(ACCELERATION_Y)) it_node->FastGetSolutionStepValue(ACCELERATION_Y) = nodal_volume_acceleration[1];
                if (!it_node->IsFixed(ACCELERATION_Z)) it_node->FastGetSolutionStepValue(ACCELERATION_Z) = nodal_volume_acceleration[2];
            }
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline CheckResult TreeContactSearch<TDim, TNumNodes>::CheckCondition(
    IndexSet::Pointer IndexesSet,
    const Condition::Pointer pCond1,
    const Condition::Pointer pCond2,
    const bool InvertedSearch
    )
{
    const IndexType index_1 = pCond1->Id();
    const IndexType index_2 = pCond2->Id();
    
    if (index_1 == index_2) // Avoiding "auto self-contact"
        return Fail;

    // Avoid conditions oriented in the same direction
    const double tolerance = 1.0e-16;
    if (norm_2(pCond1->GetValue(NORMAL) - pCond2->GetValue(NORMAL)) < tolerance)
        return Fail;

    if (pCond2->Is(SLAVE) == !InvertedSearch) { // Otherwise will not be necessary to check
        auto& indexes_set_2 = pCond2->GetValue(INDEX_SET);
        if (indexes_set_2->find(index_1) != indexes_set_2->end())
            return Fail;
    }
    
    // To avoid to repeat twice the same condition 
    if (IndexesSet->find(index_2) != IndexesSet->end())
        return AlreadyInTheMap;

    return OK;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline std::size_t TreeContactSearch<TDim, TNumNodes>::ReorderConditionsIds()
{
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
 
    std::size_t condition_id = 0;
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)  {
        ++condition_id;
        (conditions_array.begin() + i)->SetId(condition_id);
    }
    
    return condition_id;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline void TreeContactSearch<TDim, TNumNodes>::AddPotentialPairing(
    ModelPart& rComputingModelPart,
    std::size_t& ConditionId,
    Condition::Pointer pCondSlave,
    PointVector& PointsFound,
    const unsigned int NumberOfPointsFound,
    IndexSet::Pointer IndexesSet
    )
{
    // Some auxiliar values
    const double active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    const double tolerance = std::numeric_limits<double>::epsilon();
    
    // Slave geometry
    GeometryType& geom_slave = pCondSlave->GetGeometry();
    const array_1d<double, 3>& normal_slave = pCondSlave->GetValue(NORMAL);
    
    for (unsigned int i_point = 0; i_point < NumberOfPointsFound; ++i_point ) {
        bool at_least_one_node_potential_contact = false;
        
        // Master condition
        Condition::Pointer p_cond_master = PointsFound[i_point]->GetCondition();
        
        if (mCheckGap == DirectCheck) {
            // Master geometry
            const array_1d<double, 3>& normal_master = p_cond_master->GetValue(NORMAL);
            GeometryType& geom_master = p_cond_master->GetGeometry();
            
            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
                if (geom_slave[i_node].Is(ACTIVE) == false) {
                    Point projected_point;
                    double aux_distance = 0.0;
                    const array_1d<double, 3> normal = geom_slave[i_node].GetValue(NORMAL);
                    if (norm_2(normal) < tolerance)
                        aux_distance = MortarUtilities::FastProjectDirection(geom_master, geom_slave[i_node], projected_point, normal_master, normal_slave);
                    else
                        aux_distance = MortarUtilities::FastProjectDirection(geom_master, geom_slave[i_node], projected_point, normal_master, normal);
                    
                    array_1d<double, 3> result;
                    if (aux_distance <= geom_slave[i_node].FastGetSolutionStepValue(NODAL_H) * active_check_factor &&  geom_master.IsInside(projected_point, result, tolerance)) { // NOTE: This can be problematic (It depends the way IsInside() and the local_pointCoordinates() are implemented)
                        at_least_one_node_potential_contact = true;
                        geom_slave[i_node].Set(ACTIVE, true);
                    }
                    
                    aux_distance = MortarUtilities::FastProjectDirection(geom_master, geom_slave[i_node], projected_point, normal_master, -normal_master);
                    if (aux_distance <= geom_slave[i_node].FastGetSolutionStepValue(NODAL_H) * active_check_factor &&  geom_master.IsInside(projected_point, result, tolerance)) { // NOTE: This can be problematic (It depends the way IsInside() and the local_pointCoordinates() are implemented)
                        at_least_one_node_potential_contact = true;
                        geom_slave[i_node].Set(ACTIVE, true);
                    }
                }
                else
                    at_least_one_node_potential_contact = true;
            }
        } else {
            at_least_one_node_potential_contact = true;
            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                geom_slave[i_node].Set(ACTIVE, true);
        }
    
        if (at_least_one_node_potential_contact) AddPairing(rComputingModelPart, ConditionId, pCondSlave, p_cond_master, IndexesSet);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline void TreeContactSearch<TDim, TNumNodes>::CheckPairing(        
    ModelPart& rComputingModelPart,
    std::size_t& ConditionId
    )
{
    // Some auxiliar values
    const double tolerance = std::numeric_limits<double>::epsilon();
    
    // Iterate over the nodes
    ModelPart& rcontact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    NodesArrayType& nodes_array = rcontact_model_part.Nodes();
    
    // We set the auxiliar Coordinates
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        
        if (it_node->Is(MASTER) == !mInvertedSearch)
            it_node->SetValue(AUXILIAR_COORDINATES, it_node->Coordinates());
        else
            it_node->SetValue(AUXILIAR_COORDINATES, ZeroVector(3));
    }
    
    // We compute the maximal nodal h and some auxiliar values
    const double active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    // TODO: Think about this criteria
    const double distance_threshold = GetMeanNodalH(); 
//     const double distance_threshold = GetMaxNodalH(); 
//     const double distance_threshold = active_check_factor * GetMaxNodalH();
    
    // We set the mapper parameters
    Parameters mapping_parameters = Parameters(R"({"inverted_master_slave_pairing": false, "distance_threshold" : 1.0e24})" );
    mapping_parameters["inverted_master_slave_pairing"].SetBool(mThisParameters["inverted_search"].GetBool());
    mapping_parameters["distance_threshold"].SetDouble(distance_threshold);
    typedef SimpleMortarMapperProcess<TDim, TNumNodes, Variable<array_1d<double, 3>>, NonHistorical> MapperType;
    MapperType mapper = MapperType(rcontact_model_part, AUXILIAR_COORDINATES, mapping_parameters);
    mapper.Execute();
    
    // We compute now the normal gap and set the nodes under certain threshold as active
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {        
        auto it_node = nodes_array.begin() + i;

        if (it_node->Is(SLAVE) == !mInvertedSearch) {
            // We compute the gap
            const array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);
            const array_1d<double, 3>& auxiliar_coordinates = it_node->GetValue(AUXILIAR_COORDINATES);
            const array_1d<double, 3>& components_gap = ( it_node->Coordinates() - auxiliar_coordinates);
            const double gap = inner_prod(components_gap, - normal);
            it_node->SetValue(NORMAL_GAP, gap);
            
            // We activate if the node is close enough
            if (norm_2(auxiliar_coordinates) > tolerance) {
            #ifdef KRATOS_DEBUG
                KRATOS_ERROR_IF(!(it_node->SolutionStepsDataHas(NODAL_H))) << "WARNING:: NODAL_H not added" << std::endl; 
            #endif
                const double auxiliar_length = it_node->FastGetSolutionStepValue(NODAL_H) * active_check_factor;
                if (gap < auxiliar_length) it_node->Set(ACTIVE);
            }
            else
                it_node->SetValue(NORMAL_GAP, 0.0);
        }
        else
            it_node->SetValue(NORMAL_GAP, 0.0);
    }
    
    // We revert the nodes to the original position    
    if (mrMainModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY_X) == true) {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
            noalias((nodes_array.begin() + i)->Coordinates()) -= (nodes_array.begin() + i)->GetValue(DELTA_COORDINATES);
    }
    
    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(rcontact_model_part); 
    
    // Iterate in the conditions and create the new ones
    CreateAuxiliarConditions(rcontact_model_part, rComputingModelPart, ConditionId);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline void TreeContactSearch<TDim, TNumNodes>::CreateAuxiliarConditions(
    ModelPart& rContactModelPart,
    ModelPart& rComputingModelPart,
    std::size_t& ConditionId
    )
{
    // Iterate in the conditions and create the new ones
    ConditionsArrayType& conditions_array = rContactModelPart.Conditions();
    
    for(std::size_t i = 0; i < conditions_array.size(); ++i) {
        auto it_cond = conditions_array.begin() + i;
        if (it_cond->Is(SLAVE) == !mInvertedSearch) {
            IndexSet::Pointer indexes_set = it_cond->GetValue(INDEX_SET);
            for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
                Condition::Pointer p_cond_master = mrMainModelPart.pGetCondition(*it_pair); // MASTER
                AddPairing(rComputingModelPart, ConditionId, (*it_cond.base()), p_cond_master, indexes_set);
            }
        }
    } 
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
inline double TreeContactSearch<TDim, TNumNodes>::Radius(GeometryType& ThisGeometry) 
{ 
    double radius = 0.0; 
    const Point& center = ThisGeometry.Center(); 
        
    for(unsigned int i_node = 0; i_node < ThisGeometry.PointsNumber(); ++i_node)  { 
        const array_1d<double, 3>& aux_vector = center.Coordinates() - ThisGeometry[i_node].Coordinates();
        const double aux_value = inner_prod(aux_vector, aux_vector);
        if(aux_value > radius) radius = aux_value; 
    } 

    return std::sqrt(radius); 
} 

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void TreeContactSearch<TDim, TNumNodes>::ResetContactOperators()
{
    ConditionsArrayType& conditions_array = mrMainModelPart.GetSubModelPart("Contact").Conditions();
    
    #pragma omp parallel for 
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        if (it_cond->Is(SLAVE) == !mInvertedSearch)
        {            
            auto& condition_pointers = it_cond->GetValue(INDEX_SET);
            
            if (condition_pointers != nullptr)
            {
                condition_pointers->clear();
//                     condition_pointers->reserve(mAllocationSize); 
            }
        }
    }   
    
    ModelPart& computing_rcontact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact"); 
    ConditionsArrayType& computing_conditions_array = computing_rcontact_model_part.Conditions();
    const int num_computing_conditions = static_cast<int>(computing_conditions_array.size());

    #pragma omp parallel for 
    for(int i = 0; i < num_computing_conditions; ++i) 
        (computing_conditions_array.begin() + i)->Set(TO_ERASE, true);
    
    mrMainModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
SearchTreeType TreeContactSearch<TDim, TNumNodes>::ConvertSearchTree(const std::string& str)
{
    KRATOS_ERROR_IF(str == "KDOP") << "KDOP contact search: Not yet implemented" << std::endl;
    
    if(str == "InRadius") 
        return KdtreeInRadius;
    else if(str == "InBox") 
        return KdtreeInBox;
    else if (str == "KDOP")
        return Kdop;
    else
        return KdtreeInRadius;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
CheckGap TreeContactSearch<TDim, TNumNodes>::ConvertCheckGap(const std::string& str)
{
    if(str == "NoCheck") 
        return NoCheck;
    else if(str == "DirectCheck") 
        return DirectCheck;
    else if (str == "MappingCheck")
        return MappingCheck;
    else
        return MappingCheck;
}

/***********************************************************************************/
/***********************************************************************************/

template class TreeContactSearch<2, 2>;
template class TreeContactSearch<3, 3>;
template class TreeContactSearch<3, 4>;

}  // namespace Kratos.
