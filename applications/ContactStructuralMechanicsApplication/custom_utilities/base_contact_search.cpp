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
#include <algorithm>

// External includes

// Project includes
#include "input_output/logger.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/variable_utils.h"

/* Custom utilities */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/base_contact_search.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::BaseContactSearch(
    ModelPart & rMainModelPart,
    Parameters ThisParameters
    ):mrMainModelPart(rMainModelPart),
      mThisParameters(ThisParameters)
{
    KRATOS_ERROR_IF(mrMainModelPart.HasSubModelPart("Contact") == false) << "AdvancedContactSearch:: Please add the Contact submodelpart to your modelpart list" << std::endl;

    Parameters default_parameters = GetDefaultParameters();

    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    mCheckGap = this->ConvertCheckGap(mThisParameters["check_gap"].GetString());
    mInvertedSearch = mThisParameters["inverted_search"].GetBool();
    mPredefinedMasterSlave = mThisParameters["predefined_master_slave"].GetBool();

    // If we are going to consider multple searchs
    const std::string& id_name = mThisParameters["id_name"].GetString();
    mMultipleSearchs = id_name == "" ? false : true;

    // Check if the computing contact submodelpart
    const std::string sub_computing_model_part_name = "ComputingContactSub" + id_name;
    if (!(mrMainModelPart.HasSubModelPart("ComputingContact"))) { // We check if the submodelpart where the actual conditions used to compute contact are going to be computed
        ModelPart* p_computing_model_part = &mrMainModelPart.CreateSubModelPart("ComputingContact");
        p_computing_model_part->CreateSubModelPart(sub_computing_model_part_name);
    } else {
        ModelPart& r_computing_contact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact");
        if (!(r_computing_contact_model_part.HasSubModelPart(sub_computing_model_part_name)) && mMultipleSearchs) {
            r_computing_contact_model_part.CreateSubModelPart(sub_computing_model_part_name);
        } else { // We clean the existing modelpart
            ModelPart& r_sub_computing_contact_model_part = !mMultipleSearchs ? r_computing_contact_model_part : r_computing_contact_model_part.GetSubModelPart(sub_computing_model_part_name);
            ConditionsArrayType& r_conditions_array = r_sub_computing_contact_model_part.Conditions();
            VariableUtils().SetFlag(TO_ERASE, true, r_conditions_array);
            mrMainModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
        }
    }

    // Updating the base condition
    mConditionName = mThisParameters["condition_name"].GetString();
    if (mConditionName == "")
        mCreateAuxiliarConditions = false;
    else {
        mCreateAuxiliarConditions = true;
        std::ostringstream condition_name;
        condition_name << mConditionName << "Condition" << TDim << "D" << TNumNodes << "N" << mThisParameters["final_string"].GetString();
        mConditionName = condition_name.str();
    }

//     KRATOS_DEBUG_ERROR_IF_NOT(KratosComponents<Condition>::Has(mConditionName)) << "Condition " << mConditionName << " not registered" << std::endl;

    // We get the contact model part
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub" + id_name);

    // We set to zero the NORMAL_GAP
    if (mCheckGap == CheckGap::MappingCheck)
        VariableUtils().SetNonHistoricalVariable(NORMAL_GAP, 0.0, r_sub_contact_model_part.Nodes());

    // Iterate in the conditions
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();
    VariableUtils().SetFlag(ACTIVE, false, r_conditions_array);

    // We identify the type of solution
    mTypeSolution =  TypeSolution::VectorLagrangeMultiplier;
    if (mrMainModelPart.HasNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER) == false) {
        if (mrMainModelPart.HasNodalSolutionStepVariable(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)) {
            mTypeSolution = TypeSolution::NormalContactStress;
        } else {
            if (mrMainModelPart.HasNodalSolutionStepVariable(WEIGHTED_GAP)) {
                if (mrMainModelPart.Is(SLIP)) {
                    mTypeSolution = TypeSolution::FrictionalPenaltyMethod;
                } else {
                    mTypeSolution = TypeSolution::FrictionlessPenaltyMethod;
                }
            } else {
                mTypeSolution = TypeSolution::ScalarLagrangeMultiplier;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::InitializeMortarConditions()
{
    // Iterate in the conditions
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    const int num_conditions = static_cast<int>(r_conditions_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = it_cond_begin + i;
        if (!(it_cond->Has(INDEX_MAP))) {
            it_cond->SetValue(INDEX_MAP, Kratos::make_shared<IndexMap>());
//             it_cond->GetValue(INDEX_MAP)->reserve(mThisParameters["allocation_size"].GetInt());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::SetOriginDestinationModelParts(ModelPart& rModelPart)
{
    // We check if the MasterSubModelPart already exists
    const std::string& id_name = mThisParameters["id_name"].GetString();
    if (rModelPart.HasSubModelPart("MasterSubModelPart" + id_name) == false) {
        rModelPart.CreateSubModelPart("MasterSubModelPart" + id_name);
    } else {
        rModelPart.RemoveSubModelPart("MasterSubModelPart" + id_name);
        rModelPart.CreateSubModelPart("MasterSubModelPart" + id_name);
    }
    // We check if the SlaveSubModelPart already exists
    if (rModelPart.HasSubModelPart("SlaveSubModelPart" + id_name) == false) {
        rModelPart.CreateSubModelPart("SlaveSubModelPart" + id_name);
    } else {
        rModelPart.RemoveSubModelPart("SlaveSubModelPart" + id_name);
        rModelPart.CreateSubModelPart("SlaveSubModelPart" + id_name);
    }

    ModelPart& r_master_model_part = rModelPart.GetSubModelPart("MasterSubModelPart" + id_name);
    ModelPart& r_slave_model_part = rModelPart.GetSubModelPart("SlaveSubModelPart" + id_name);

    // The vectors containing the ids
    std::vector<IndexType> slave_nodes_ids,  master_nodes_ids;
    std::vector<IndexType> slave_conditions_ids, master_conditions_ids;

    // Begin iterators
    const auto it_node_begin = rModelPart.NodesBegin();
    const auto it_cond_begin = rModelPart.ConditionsBegin();

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        std::vector<IndexType> slave_nodes_ids_buffer, master_nodes_ids_buffer;
        std::vector<IndexType> slave_conditions_ids_buffer, master_conditions_ids_buffer;

        #pragma omp for
        for(int i=0; i<static_cast<int>(rModelPart.Nodes().size()); ++i) {
            auto it_node = it_node_begin + i;

            if (it_node->Is(SLAVE) == !mInvertedSearch) {
                slave_nodes_ids_buffer.push_back(it_node->Id());
            }
            if (it_node->Is(MASTER) == !mInvertedSearch) {
                master_nodes_ids_buffer.push_back(it_node->Id());
            }
        }

        #pragma omp for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;

            if (it_cond->Is(SLAVE) == !mInvertedSearch) {
                slave_conditions_ids_buffer.push_back(it_cond->Id());
            }
            if (it_cond->Is(MASTER) == !mInvertedSearch) {
                master_conditions_ids_buffer.push_back(it_cond->Id());
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(slave_nodes_ids_buffer.begin(),slave_nodes_ids_buffer.end(),back_inserter(slave_nodes_ids));
            std::move(master_nodes_ids_buffer.begin(),master_nodes_ids_buffer.end(),back_inserter(master_nodes_ids));
            std::move(slave_conditions_ids_buffer.begin(),slave_conditions_ids_buffer.end(),back_inserter(slave_conditions_ids));
            std::move(master_conditions_ids_buffer.begin(),master_conditions_ids_buffer.end(),back_inserter(master_conditions_ids));
        }
    }

    // Finally we add the nodes and conditions to the submodelparts
    r_slave_model_part.AddNodes(slave_nodes_ids);
    r_slave_model_part.AddConditions(slave_conditions_ids);
    r_master_model_part.AddNodes(master_nodes_ids);
    r_master_model_part.AddConditions(master_conditions_ids);

    KRATOS_ERROR_IF(r_master_model_part.Conditions().size() == 0) << "No origin conditions. Check your flags are properly set" << std::endl;
    KRATOS_ERROR_IF(r_slave_model_part.Conditions().size() == 0) << "No destination conditions. Check your flags are properly set" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ClearMortarConditions()
{
    ResetContactOperators();

    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();

    switch(mTypeSolution) {
        case TypeSolution::VectorLagrangeMultiplier :
            ClearComponentsMortarConditions(r_nodes_array);
            break;
        case TypeSolution::ScalarLagrangeMultiplier :
            ClearScalarMortarConditions(r_nodes_array);
            break;
        case TypeSolution::NormalContactStress :
            ClearALMFrictionlessMortarConditions(r_nodes_array);
            break;
        case TypeSolution::FrictionlessPenaltyMethod :
            break;
        case TypeSolution::FrictionalPenaltyMethod :
            break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckContactModelParts()
{
    // Iterate in the conditions
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();

    const SizeType total_number_conditions = mrMainModelPart.GetRootModelPart().NumberOfConditions();

    std::vector<Condition::Pointer> auxiliar_conditions_vector;

    for(Condition& r_cond : r_conditions_array) {
        if (r_cond.Is(MARKER)) {
            // Setting the flag to remove
            r_cond.Set(TO_ERASE, true);

            // Creating new condition
            Condition::Pointer p_new_cond = r_cond.Clone(total_number_conditions + r_cond.Id(), r_cond.GetGeometry());
            auxiliar_conditions_vector.push_back(p_new_cond);

            p_new_cond->SetData(r_cond.GetData()); // TODO: Remove when fixed on the core
            p_new_cond->SetValue(INDEX_MAP, Kratos::make_shared<IndexMap>());
//             p_new_cond->GetValue(INDEX_MAP)->clear();
//             p_new_cond->GetValue(INDEX_MAP)->reserve(mThisParameters["allocation_size"].GetInt());
            p_new_cond->Set(Flags(r_cond));
            p_new_cond->Set(MARKER, true);
        } else {
            // Setting the flag to mark
            r_cond.Set(MARKER, true);
        }
    }

    // Finally we add the new conditions to the model part
    r_sub_contact_model_part.RemoveConditions(TO_ERASE);
    // Reorder ids (in order to keep the ids consistent)
    for (int i = 0; i < static_cast<int>(auxiliar_conditions_vector.size()); ++i) {
        auxiliar_conditions_vector[i]->SetId(total_number_conditions + i + 1);
    }
    ConditionsArrayType aux_conds;
    aux_conds.GetContainer() = auxiliar_conditions_vector;
    r_sub_contact_model_part.AddConditions(aux_conds.begin(), aux_conds.end());

    // Unsetting TO_ERASE
    VariableUtils().SetFlag(TO_ERASE, false, r_contact_model_part.Conditions());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CreatePointListMortar()
{
    // Clearing the vector
    mPointListDestination.clear();

    // Iterate in the conditions
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    for(IndexType i = 0; i < r_conditions_array.size(); ++i) {
        auto it_cond = it_cond_begin + i;
        if (it_cond->Is(MASTER) == !mInvertedSearch || !mPredefinedMasterSlave) {
            mPointListDestination.push_back(Kratos::make_shared<PointItem>((*it_cond.base())));
        }
    }

#ifdef KRATOS_DEBUG
    // NOTE: We check the list
    for (IndexType i_point = 0; i_point < mPointListDestination.size(); ++i_point )
        mPointListDestination[i_point]->Check();
#endif
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::UpdatePointListMortar()
{
    // We check if we are in a dynamic or static case
    const bool dynamic = mThisParameters["dynamic_search"].GetBool() ? mrMainModelPart.HasNodalSolutionStepVariable(VELOCITY) : false;
    const double delta_time = (dynamic) ? mrMainModelPart.GetProcessInfo()[DELTA_TIME] : 0.0;

    // The contact model parts
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());

    // We compute the delta displacement
    if (dynamic) {
        ContactUtilities::ComputeStepJump(r_sub_contact_model_part, delta_time);
    }

    if (mCheckGap == CheckGap::MappingCheck && dynamic) {
        NodesArrayType& r_update_r_nodes_array = r_sub_contact_model_part.Nodes();
        const auto it_node_begin = r_update_r_nodes_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_update_r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
            noalias(it_node->Coordinates()) += it_node->GetValue(DELTA_COORDINATES);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mPointListDestination.size()); ++i)
        mPointListDestination[i]->UpdatePoint();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::UpdateMortarConditions()
{
    // We update the list of points
    UpdatePointListMortar();

    // The contact model parts
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(r_sub_contact_model_part);

    // We get the computing model part
    IndexType condition_id = GetMaximumConditionsIds();
    const std::string sub_computing_model_part_name = "ComputingContactSub" + mThisParameters["id_name"].GetString();
    ModelPart& r_computing_contact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact");
    ModelPart& r_sub_computing_contact_model_part = !mMultipleSearchs ? r_computing_contact_model_part : r_computing_contact_model_part.GetSubModelPart(sub_computing_model_part_name);

    // We reset the computing contact model part in case of already initialized
    if (r_sub_computing_contact_model_part.Conditions().size() > 0)
        ClearMortarConditions();

    // We check if we are in a dynamic or static case
    const bool dynamic = mThisParameters["dynamic_search"].GetBool() ? mrMainModelPart.HasNodalSolutionStepVariable(VELOCITY) : false;

    // Some auxiliar values
    const IndexType allocation_size = mThisParameters["allocation_size"].GetInt();              // Allocation size for the vectors and max number of potential results
    const double search_factor = mThisParameters["search_factor"].GetDouble();                  // The search factor to be considered
    SearchTreeType type_search = ConvertSearchTree(mThisParameters["type_search"].GetString()); // The search tree considered
    IndexType bucket_size = mThisParameters["bucket_size"].GetInt();                            // Bucket size for kd-tree

    // Create a tree
    // It will use a copy of mNodeList (a std::vector which contains pointers)
    // Copying the list is required because the tree will reorder it for efficiency
    KDTree tree_points(mPointListDestination.begin(), mPointListDestination.end(), bucket_size);

    // Auxiliar model parts and components
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();
    const int num_conditions = static_cast<int>(r_conditions_array.size());

    // In case of not predefined master/slave we reset the flags
    if (!mPredefinedMasterSlave) {
        NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
        VariableUtils().SetFlag(SLAVE, false, r_nodes_array);
        VariableUtils().SetFlag(MASTER, false, r_nodes_array);
        VariableUtils().SetFlag(SLAVE, false, r_conditions_array);
        VariableUtils().SetFlag(MASTER, false, r_conditions_array);
    }

    // Now we iterate over the conditions
//     #pragma omp parallel for firstprivate(tree_points) // FIXME: Make me parallel!!!
    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = r_conditions_array.begin() + i;

        if (!mPredefinedMasterSlave || it_cond->Is(SLAVE) == !mInvertedSearch) {
            // Initialize values
            PointVector points_found(allocation_size);

            IndexType number_points_found = 0;

            if (type_search == SearchTreeType::KdtreeInRadius) {
                GeometryType& r_geometry = it_cond->GetGeometry();
                const Point& r_center = dynamic ? Point(ContactUtilities::GetHalfJumpCenter(r_geometry)) : r_geometry.Center(); // NOTE: Center in half delta time or real center

                const double search_radius = search_factor * Radius(it_cond->GetGeometry());

                number_points_found = tree_points.SearchInRadius(r_center, search_radius, points_found.begin(), allocation_size);
            } else if (type_search == SearchTreeType::KdtreeInBox) {
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

                const array_1d<double,3> normal_min = MortarUtilities::GaussPointUnitNormal(N_min, it_cond->GetGeometry());
                const array_1d<double,3> normal_max = MortarUtilities::GaussPointUnitNormal(N_max, it_cond->GetGeometry());

                ContactUtilities::ScaleNode<NodeType>(min_point, normal_min, length_search);
                ContactUtilities::ScaleNode<NodeType>(max_point, normal_max, length_search);

                number_points_found = tree_points.SearchInBox(min_point, max_point, points_found.begin(), allocation_size);
            } else
                KRATOS_ERROR << " The type search is not implemented yet does not exist!!!!. SearchTreeType = " << mThisParameters["type_search"].GetString() << std::endl;

            if (number_points_found > 0) {
                // We resize the vector to the actual size
//                 points_found.resize(number_points_found); // NOTE: May be ineficient

            #ifdef KRATOS_DEBUG
                // NOTE: We check the list
                for (IndexType i_point = 0; i_point < number_points_found; ++i_point )
                    points_found[i_point]->Check();
//                 KRATOS_INFO("Check search") << "The search is properly done" << std::endl;
            #endif

                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue(INDEX_MAP);

                // If not active we check if can be potentially in contact
                if (mCheckGap == CheckGap::MappingCheck) {
                    for (IndexType i_point = 0; i_point < number_points_found; ++i_point ) {
                        Condition::Pointer p_cond_master = points_found[i_point]->GetCondition();
                        const CheckResult condition_checked_right = CheckCondition(p_indexes_pairs, (*it_cond.base()), p_cond_master, mInvertedSearch);

                        if (condition_checked_right == CheckResult::OK)
                            p_indexes_pairs->AddId(p_cond_master->Id());
                    }
                } else {
                    AddPotentialPairing(r_sub_computing_contact_model_part, condition_id, (*it_cond.base()), points_found, number_points_found, p_indexes_pairs);
                }
            }
        }
    }

    // In case of not predefined master/slave we assign the master/slave nodes and conditions
    if (!mPredefinedMasterSlave)
        NotPredefinedMasterSlave(r_sub_contact_model_part);

    // We create the submodelparts for master and slave
    SetOriginDestinationModelParts(r_sub_contact_model_part);

    // We map the Coordinates to the slave side from the master
    if (mCheckGap == CheckGap::MappingCheck) {
        CheckPairing(r_sub_computing_contact_model_part, condition_id);
    } else {
        // We revert the nodes to the original position
        if (mThisParameters["dynamic_search"].GetBool()) {
            if (mrMainModelPart.HasNodalSolutionStepVariable(VELOCITY)) {
                NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
                #pragma omp parallel for
                for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                    auto it_node = r_nodes_array.begin() + i;
                    noalias(it_node->Coordinates()) -= it_node->GetValue(DELTA_COORDINATES);
                }
            }
        }
        // We compute the weighted reaction
        ComputeWeightedReaction();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::AddPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId,
    Condition::Pointer pCondSlave,
    Condition::Pointer pCondMaster
    )
{
    if (mCreateAuxiliarConditions) { // We add the ID and we create a new auxiliar condition
        ++rConditionId;
        Condition::Pointer p_auxiliar_condition = rComputingModelPart.CreateNewCondition(mConditionName, rConditionId, pCondSlave->GetGeometry(), pCondSlave->pGetProperties());
        // We set the geometrical values
        IndexMap::Pointer ids_destination = pCondSlave->GetValue(INDEX_MAP);
        ids_destination->SetNewEntityId(pCondMaster->Id(), rConditionId);
        p_auxiliar_condition->SetValue(PAIRED_GEOMETRY, pCondMaster->pGetGeometry());
        p_auxiliar_condition->SetValue(NORMAL, pCondSlave->GetValue(NORMAL));
        p_auxiliar_condition->SetValue(PAIRED_NORMAL, pCondMaster->GetValue(NORMAL));
        // We activate the condition and initialize it
        p_auxiliar_condition->Set(ACTIVE, true);
        p_auxiliar_condition->Initialize();
        // TODO: Check this!!
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::AddPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId,
    Condition::Pointer pCondSlave,
    Condition::Pointer pCondMaster,
    IndexMap::Pointer IndexesPairs
    )
{
    IndexesPairs->AddId(pCondMaster->Id());

    AddPairing(rComputingModelPart, rConditionId, pCondSlave, pCondMaster);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckMortarConditions()
{
    // Iterate in the conditions
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();

    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = r_conditions_array.begin() + i;

        if (it_cond->Has(INDEX_MAP)) {
            IndexMap::Pointer ids_destination = it_cond->GetValue(INDEX_MAP);
            if (ids_destination->size() > 0) {
                KRATOS_INFO("Check paired conditions (Origin)") << "Origin condition ID:" << it_cond->Id() << " Number of pairs: " << ids_destination->size() << std::endl;
                KRATOS_INFO("Check paired conditions (Destination)") << ids_destination->Info();
            }
        }
    }

    // Iterate over the nodes
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        if (it_node->Is(ACTIVE))
            KRATOS_INFO("Check paired nodes") << "Node: " << it_node->Id() << " is active" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::InvertSearch()
{
    mInvertedSearch = !mInvertedSearch;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ClearScalarMortarConditions(NodesArrayType& rNodesArray)
{
    VariableUtils().SetVariableForFlag(SCALAR_LAGRANGE_MULTIPLIER, 0.0, rNodesArray, ACTIVE, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ClearComponentsMortarConditions(NodesArrayType& rNodesArray)
{
    const array_1d<double, 3> zero_array = ZeroVector(3);
    VariableUtils().SetVariableForFlag(VECTOR_LAGRANGE_MULTIPLIER, zero_array, rNodesArray, ACTIVE, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ClearALMFrictionlessMortarConditions(NodesArrayType& rNodesArray)
{
    VariableUtils().SetVariableForFlag(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0.0, rNodesArray, ACTIVE, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline typename BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckResult BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckCondition(
    IndexMap::Pointer pIndexesPairs,
    const Condition::Pointer pCond1,
    const Condition::Pointer pCond2,
    const bool InvertedSearch
    )
{
    const IndexType index_1 = pCond1->Id();
    const IndexType index_2 = pCond2->Id();

    if (index_1 == index_2) // Avoiding "auto self-contact"
        return CheckResult::Fail;

    // Avoid conditions oriented in the same direction
    const double tolerance = 1.0e-16;
    if (norm_2(pCond1->GetValue(NORMAL) - pCond2->GetValue(NORMAL)) < tolerance)
        return CheckResult::Fail;

    // Otherwise will not be necessary to check
    if (!mPredefinedMasterSlave || pCond2->Is(SLAVE) == !InvertedSearch) {
        auto p_indexes_pairs_2 = pCond2->GetValue(INDEX_MAP);
        if (p_indexes_pairs_2->find(index_1) != p_indexes_pairs_2->end())
            return CheckResult::Fail;
    }

    // To avoid to repeat twice the same condition
    if (pIndexesPairs->find(index_2) != pIndexesPairs->end())
        return CheckResult::AlreadyInTheMap;

    return CheckResult::OK;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::NotPredefinedMasterSlave(ModelPart& rModelPart)
{
    // We iterate over the conditions
    ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    const int num_conditions = static_cast<int>(r_conditions_array.size());

    std::vector<IndexType> master_conditions_ids;

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        std::vector<IndexType> master_conditions_ids_buffer;

        #pragma omp for
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = it_cond_begin + i;
            IndexMap::Pointer p_indexes_pairs = it_cond->GetValue(INDEX_MAP);
            if (p_indexes_pairs->size() > 0) {
                it_cond->Set(SLAVE, true);
                for (auto& i_pair : *p_indexes_pairs) {
                    master_conditions_ids_buffer.push_back(i_pair.first);
                }
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(master_conditions_ids_buffer.begin(),master_conditions_ids_buffer.end(),back_inserter(master_conditions_ids));
        }
    }

    // We create an auxiliar model part to add the MASTER flag
    rModelPart.CreateSubModelPart("AuxMasterModelPart");
    ModelPart& aux_model_part = rModelPart.GetSubModelPart("AuxMasterModelPart");

    // Remove duplicates
    std::sort( master_conditions_ids.begin(), master_conditions_ids.end() );
    master_conditions_ids.erase( std::unique( master_conditions_ids.begin(), master_conditions_ids.end() ), master_conditions_ids.end() );

    // Add to the auxiliar model part
    aux_model_part.AddConditions(master_conditions_ids);

    // Set the flag
    VariableUtils().SetFlag(MASTER, true, aux_model_part.Conditions());

    // Remove auxiliar model part
    rModelPart.RemoveSubModelPart("AuxMasterModelPart");

    // Now we iterate over the conditions to set the nodes indexes
    #pragma omp parallel for
    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = r_conditions_array.begin() + i;
        if (it_cond->Is(SLAVE)) {
            GeometryType& r_geometry = it_cond->GetGeometry();
            for (NodeType& r_node : r_geometry) {
                r_node.SetLock();
                r_node.Set(SLAVE, true);
                r_node.UnSetLock();
            }
        }
        if (it_cond->Is(MASTER)) {
            GeometryType& r_geometry = it_cond->GetGeometry();
            for (NodeType& r_node : r_geometry) {
                r_node.SetLock();
                r_node.Set(MASTER, true);
                r_node.UnSetLock();
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline IndexType BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::GetMaximumConditionsIds()
{
    ConditionsArrayType& r_conditions_array = mrMainModelPart.Conditions();

    IndexType condition_id = 0;
    for(IndexType i = 0; i < r_conditions_array.size(); ++i)  {
        auto it_cond = r_conditions_array.begin() + i;
        const IndexType id = it_cond->GetId();
        if (id > condition_id)
            condition_id = id;
    }

    return condition_id;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::AddPotentialPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId,
    Condition::Pointer pCondSlave,
    PointVector& rPointsFound,
    const IndexType NumberOfPointsFound,
    IndexMap::Pointer IndexesPairs
    )
{
    // Some auxiliar values
    const double active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    const bool frictional_problem = mrMainModelPart.Is(SLIP);

    // Slave geometry
    GeometryType& r_geom_slave = pCondSlave->GetGeometry();
    const array_1d<double, 3>& r_normal_slave = pCondSlave->GetValue(NORMAL);

    for (IndexType i_point = 0; i_point < NumberOfPointsFound; ++i_point ) {
        bool at_least_one_node_potential_contact = false;

        // Master condition
        Condition::Pointer p_cond_master = rPointsFound[i_point]->GetCondition();

        if (mCheckGap == CheckGap::DirectCheck) {
            // Master geometry
            const array_1d<double, 3>& normal_master = p_cond_master->GetValue(NORMAL);
            GeometryType& r_geom_master = p_cond_master->GetGeometry();

            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                if (r_geom_slave[i_node].IsNot(ACTIVE)) {
                    Point projected_point;
                    double aux_distance = 0.0;
                    const array_1d<double, 3> normal = r_geom_slave[i_node].GetValue(NORMAL);
                    if (norm_2(normal) < ZeroTolerance)
                        aux_distance = GeometricalProjectionUtilities::FastProjectDirection(r_geom_master, r_geom_slave[i_node], projected_point, normal_master, r_normal_slave);
                    else
                        aux_distance = GeometricalProjectionUtilities::FastProjectDirection(r_geom_master, r_geom_slave[i_node], projected_point, normal_master, normal);

                    array_1d<double, 3> result;
                    if (aux_distance <= r_geom_slave[i_node].FastGetSolutionStepValue(NODAL_H) * active_check_factor &&  r_geom_master.IsInside(projected_point, result, ZeroTolerance)) { // NOTE: This can be problematic (It depends the way IsInside() and the local_pointCoordinates() are implemented)
                        at_least_one_node_potential_contact = true;
                        r_geom_slave[i_node].Set(ACTIVE, true);
                        if (mTypeSolution == TypeSolution::VectorLagrangeMultiplier && frictional_problem) {
                            NodeType& node = r_geom_slave[i_node];
                            if (norm_2(node.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
                                if (node.GetValue(FRICTION_COEFFICIENT) < ZeroTolerance) {
                                    node.Set(SLIP, true);
                                } else {
                                    node.Set(SLIP, false);
                                }
                            }
                        }
                    }

                    aux_distance = GeometricalProjectionUtilities::FastProjectDirection(r_geom_master, r_geom_slave[i_node], projected_point, normal_master, -normal_master);
                    if (aux_distance <= r_geom_slave[i_node].FastGetSolutionStepValue(NODAL_H) * active_check_factor &&  r_geom_master.IsInside(projected_point, result, ZeroTolerance)) { // NOTE: This can be problematic (It depends the way IsInside() and the local_pointCoordinates() are implemented)
                        at_least_one_node_potential_contact = true;
                        r_geom_slave[i_node].Set(ACTIVE, true);
                        if (mTypeSolution == TypeSolution::VectorLagrangeMultiplier && frictional_problem) {
                            NodeType& node = r_geom_slave[i_node];
                            if (norm_2(node.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
                                if (node.GetValue(FRICTION_COEFFICIENT) < ZeroTolerance) {
                                    node.Set(SLIP, true);
                                } else {
                                    node.Set(SLIP, false);
                                }
                            }
                        }
                    }
                } else
                    at_least_one_node_potential_contact = true;
            }
        } else {
            at_least_one_node_potential_contact = true;
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                r_geom_slave[i_node].Set(ACTIVE, true);
                if (mTypeSolution == TypeSolution::VectorLagrangeMultiplier && frictional_problem) {
                    NodeType& r_node = r_geom_slave[i_node];
                    if (norm_2(r_geom_slave[i_node].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
                        if (norm_2(r_node.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
                            if (r_node.GetValue(FRICTION_COEFFICIENT) < ZeroTolerance) {
                                r_node.Set(SLIP, true);
                            } else {
                                r_node.Set(SLIP, false);
                            }
                        }
                    }
                }
            }
        }

        if (at_least_one_node_potential_contact)
            AddPairing(rComputingModelPart, rConditionId, pCondSlave, p_cond_master, IndexesPairs);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId
    )
{
    // Getting the corresponding submodelparts
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());

    // We set the gap to an enormous value in order to initialize it
    VariableUtils().SetNonHistoricalVariable(NORMAL_GAP, 1.0e12, r_sub_contact_model_part.Nodes());

    // We compute the gap in the slave
    ComputeMappedGap(!mInvertedSearch);

    // We revert the nodes to the original position
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
    if (mThisParameters["dynamic_search"].GetBool()) {
        if (mrMainModelPart.HasNodalSolutionStepVariable(VELOCITY)) {
            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = r_nodes_array.begin() + i;
                noalias(it_node->Coordinates()) -= it_node->GetValue(DELTA_COORDINATES);
            }
        }
    }

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(r_sub_contact_model_part);

    // Iterate in the conditions and create the new ones
    CreateAuxiliarConditions(r_sub_contact_model_part, rComputingModelPart, rConditionId);

    // We compute the weighted reaction
    ComputeWeightedReaction();

    // Finally we compute the active/inactive nodes
    ComputeActiveInactiveNodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ComputeMappedGap(const bool SearchOrientation)
{
    // We get the process info
    const ProcessInfo& r_process_info = mrMainModelPart.GetProcessInfo();

    // Iterate over the nodes
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ModelPart& r_master_model_part = r_sub_contact_model_part.GetSubModelPart("MasterSubModelPart"+mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array_master = r_master_model_part.Nodes();
    const auto it_node_begin_master = r_nodes_array_master.begin();
    ModelPart& r_slave_model_part = r_sub_contact_model_part.GetSubModelPart("SlaveSubModelPart"+mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array_slave = r_slave_model_part.Nodes();
    const auto it_node_begin_slave = r_nodes_array_slave.begin();

    // We set the auxiliar Coordinates
    const array_1d<double, 3> zero_array = ZeroVector(3);
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array_master.size()); ++i) {
        auto it_node = it_node_begin_master + i;

        if (SearchOrientation) {
            it_node->SetValue(AUXILIAR_COORDINATES, it_node->Coordinates());
        } else {
            it_node->SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array_slave.size()); ++i) {
        auto it_node = it_node_begin_slave + i;

        if (!SearchOrientation) {
            it_node->SetValue(AUXILIAR_COORDINATES, it_node->Coordinates());
        } else {
            it_node->SetValue(AUXILIAR_COORDINATES, zero_array);
        }
    }

    // Switch MASTER/SLAVE
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
    if (!SearchOrientation)
        SwitchFlagNodes(r_nodes_array);

    // We set the mapper parameters
    Parameters mapping_parameters = Parameters(R"({"distance_threshold" : 1.0e24, "remove_isolated_conditions" : true, "origin_variable_historical" : false, "destination_variable_historical" : false})" );
    if (r_process_info.Has(DISTANCE_THRESHOLD)) {
        mapping_parameters["distance_threshold"].SetDouble(r_process_info[DISTANCE_THRESHOLD]);
    }
    MapperType mapper(r_master_model_part, r_slave_model_part, AUXILIAR_COORDINATES, mapping_parameters);
    mapper.Execute();

    // Switch again MASTER/SLAVE
    if (!SearchOrientation)
        SwitchFlagNodes(r_nodes_array);

    // We compute now the normal gap and set the nodes under certain threshold as active
    array_1d<double, 3> normal, auxiliar_coordinates, components_gap;
    double gap;
    const auto it_node_begin = r_nodes_array.begin();
    #pragma omp parallel for firstprivate(gap, normal, auxiliar_coordinates, components_gap)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        if (it_node->Is(SLAVE) == SearchOrientation) {
            // We compute the gap
            noalias(normal) = it_node->FastGetSolutionStepValue(NORMAL);
            noalias(auxiliar_coordinates) = it_node->GetValue(AUXILIAR_COORDINATES);
            noalias(components_gap) = ( it_node->Coordinates() - auxiliar_coordinates);
            gap = inner_prod(components_gap, - normal);

            // We activate if the node is close enough
            if (norm_2(auxiliar_coordinates) > ZeroTolerance)
                it_node->SetValue(NORMAL_GAP, gap);
        } else {
            it_node->SetValue(NORMAL_GAP, 0.0);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ComputeActiveInactiveNodes()
{
    // We get the process info
    const ProcessInfo& r_process_info = mrMainModelPart.GetProcessInfo();

    // The penalty value and scale factor
    const double common_epsilon = r_process_info[INITIAL_PENALTY];
    const double scale_factor = r_process_info[SCALE_FACTOR];

    // Iterate over the nodes
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();

    // We compute now the normal gap and set the nodes under certain threshold as active
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = r_nodes_array.begin() + i;
        if (it_node->Is(SLAVE) == !mInvertedSearch) {
//             if (it_node->GetValue(NORMAL_GAP) < ZeroTolerance) {
            if (it_node->GetValue(NORMAL_GAP) < GapThreshold * it_node->FastGetSolutionStepValue(NODAL_H)) {
                SetActiveNode(it_node, common_epsilon, scale_factor);
            } else {
            #ifdef KRATOS_DEBUG
                KRATOS_WARNING_IF("BaseContactSearch", it_node->Is(ACTIVE)) << "WARNING: A node that used to be active is not active anymore. Check that. Node ID: " << it_node->Id() << std::endl;
            #endif
                SetInactiveNode(it_node);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::SetActiveNode(
    NodesArrayType::iterator ItNode,
    const double CommonEpsilon,
    const double ScaleFactor
    )
{
    // We activate
    ItNode->Set(ACTIVE, true);
    ItNode->Set(MARKER, true);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::SetInactiveNode(NodesArrayType::iterator ItNode)
{
    // If the node has been already actived we do not inactivate
    if (ItNode->IsNot(MARKER)) {
        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        if (ItNode->Is(ACTIVE) ) {
            ItNode->Set(ACTIVE, false);
            switch(mTypeSolution) {
                case TypeSolution::VectorLagrangeMultiplier :
                    noalias(ItNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) = zero_array;
                    break;
                case TypeSolution::ScalarLagrangeMultiplier :
                    ItNode->FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER) = 0.0;
                    break;
                case TypeSolution::NormalContactStress :
                    ItNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) = 0.0;
                    break;
                case TypeSolution::FrictionlessPenaltyMethod :
                    break;
                case TypeSolution::FrictionalPenaltyMethod :
                    break;
            }
        }

        // We set the gap to zero (in order to have something "visible" to post process)
        ItNode->SetValue(NORMAL_GAP, 0.0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ComputeWeightedReaction()
{
    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // Auxiliar gap
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();
    switch(mTypeSolution) {
        case TypeSolution::VectorLagrangeMultiplier :
            if (mrMainModelPart.Is(SLIP)) {
                VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
                VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_array, r_nodes_array);
            } else if (mrMainModelPart.Is(CONTACT)) {
                VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
            } else
                VariableUtils().SetVectorVar(WEIGHTED_VECTOR_RESIDUAL, zero_array, r_nodes_array);
            break;
        case TypeSolution::ScalarLagrangeMultiplier :
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_SCALAR_RESIDUAL, 0.0, r_nodes_array);
            break;
        case TypeSolution::NormalContactStress :
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
            break;
        case TypeSolution::FrictionlessPenaltyMethod :
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
            break;
        case TypeSolution::FrictionalPenaltyMethod :
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_nodes_array);
            VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_array, r_nodes_array);
            break;
    }

    // Compute explicit contibution of the conditions
    const std::string sub_computing_model_part_name = "ComputingContactSub" + mThisParameters["id_name"].GetString();
    ModelPart& r_computing_contact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact");
    ModelPart& r_sub_computing_contact_model_part = !mMultipleSearchs ? r_computing_contact_model_part : r_computing_contact_model_part.GetSubModelPart(sub_computing_model_part_name);
    ContactUtilities::ComputeExplicitContributionConditions(r_sub_computing_contact_model_part);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CreateAuxiliarConditions(
    ModelPart& rContactModelPart,
    ModelPart& rComputingModelPart,
    IndexType& rConditionId
    )
{
    // Iterate in the conditions and create the new ones
    ConditionsArrayType& r_conditions_array = rContactModelPart.Conditions();

    // In case of debug mode
    if (mThisParameters["debug_mode"].GetBool()) {
        std::filebuf debug_buffer;
        debug_buffer.open("original_conditions_normal_debug_" + rContactModelPart.Name() + "_step=" + std::to_string( rContactModelPart.GetProcessInfo()[STEP]) + ".out",std::ios::out);
        std::ostream os(&debug_buffer);
        for (const auto& cond : r_conditions_array) {
            const array_1d<double, 3>& r_normal = cond.GetValue(NORMAL);
            os << "Condition " << cond.Id() << "\tNodes ID:";
            for (auto& r_node : cond.GetGeometry()) {
                os << "\t" << r_node.Id();
            }
            os << "\tNORMAL: " << r_normal[0] << "\t" << r_normal[1] << "\t" << r_normal[2] <<"\n";
        }
        debug_buffer.close();
    }

    // Actually creating the new conditions
    for(IndexType i = 0; i < r_conditions_array.size(); ++i) {
        auto it_cond = r_conditions_array.begin() + i;
        if (it_cond->Is(SLAVE) == !mInvertedSearch) {
            IndexMap::Pointer p_indexes_pairs = it_cond->GetValue(INDEX_MAP);
            for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                if (it_pair->second == 0) { // If different than 0 it is an existing condition
                    Condition::Pointer p_cond_master = mrMainModelPart.pGetCondition(it_pair->first); // MASTER
                    AddPairing(rComputingModelPart, rConditionId, (*it_cond.base()), p_cond_master);
                }
            }
        }
    }

    // In case of debug mode
    if (mThisParameters["debug_mode"].GetBool()) {
        std::filebuf debug_buffer;
        debug_buffer.open("created_conditions_normal_debug_" + rContactModelPart.Name() + "_step=" + std::to_string( rContactModelPart.GetProcessInfo()[STEP]) + ".out",std::ios::out);
        std::ostream os(&debug_buffer);
        for (const auto& cond : rComputingModelPart.Conditions()) {
            const array_1d<double, 3>& r_normal = cond.GetValue(NORMAL);
            os << "Condition " << cond.Id() << "\tNodes ID:";
            for (auto& r_node : cond.GetGeometry()) {
                os << "\t" << r_node.Id();
            }
            os << "\tNORMAL: " << r_normal[0] << "\t" << r_normal[1] << "\t" << r_normal[2] <<"\n";
        }
        debug_buffer.close();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline double BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::Radius(GeometryType& ThisGeometry)
{
    double radius = 0.0;
    const Point& r_center = ThisGeometry.Center();

    for(IndexType i_node = 0; i_node < ThisGeometry.PointsNumber(); ++i_node)  {
        const array_1d<double, 3>& aux_vector = r_center.Coordinates() - ThisGeometry[i_node].Coordinates();
        const double aux_value = norm_2(aux_vector);
        if(aux_value > radius)
            radius = aux_value;
    }

    return radius;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ResetContactOperators()
{
    // We iterate over the conditions
    ModelPart& r_contact_model_part = mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = !mMultipleSearchs ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+mThisParameters["id_name"].GetString());
    ConditionsArrayType& r_conditions_array = r_sub_contact_model_part.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    if (mrMainModelPart.Is(MODIFIED)) { // It has been remeshed. We remove everything

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            if (it_cond->Is(SLAVE) == !mInvertedSearch) {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue(INDEX_MAP);

                if (p_indexes_pairs != nullptr) {
                    p_indexes_pairs->clear();
//                     p_indexes_pairs->reserve(mAllocationSize);
                }
            }
        }

        // We remove all the computing conditions conditions
        const std::string sub_computing_model_part_name = "ComputingContactSub" + mThisParameters["id_name"].GetString();
        ModelPart& r_computing_contact_model_part = mrMainModelPart.GetSubModelPart("ComputingContact");
        ModelPart& r_sub_computing_contact_model_part = !mMultipleSearchs ? r_computing_contact_model_part : r_computing_contact_model_part.GetSubModelPart(sub_computing_model_part_name);
        ConditionsArrayType& r_computing_conditions_array = r_sub_computing_contact_model_part.Conditions();
        VariableUtils().SetFlag(TO_ERASE, true, r_computing_conditions_array);
    } else {
        // We iterate, but not in OMP
        for(IndexType i = 0; i < r_conditions_array.size(); ++i) {
            auto it_cond = r_conditions_array.begin() + i;
            if (it_cond->Is(SLAVE) == !mInvertedSearch) {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue(INDEX_MAP);
                if (p_indexes_pairs != nullptr) {
                    // The vector with the ids to remove
                    std::vector<IndexType> inactive_conditions_ids;
                    for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                        Condition::Pointer p_cond = mrMainModelPart.pGetCondition(it_pair->second);
                        if (p_cond->IsNot(ACTIVE)) {
                            p_cond->Set(TO_ERASE, true);
                            inactive_conditions_ids.push_back(it_pair->first);
                        }
                    }
                    for (auto& i_to_remove : inactive_conditions_ids) {
                        p_indexes_pairs->RemoveId(i_to_remove);
                    }
                }
            }
        }
    }

    mrMainModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
typename BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::SearchTreeType BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ConvertSearchTree(const std::string& str)
{
    KRATOS_ERROR_IF(str == "KDOP") << "KDOP contact search: Not yet implemented" << std::endl;

    if(str == "InRadius" || str == "in_radius")
        return SearchTreeType::KdtreeInRadius;
    else if(str == "InBox" || str == "in_box")
        return SearchTreeType::KdtreeInBox;
    else if (str == "KDOP" || str == "kdop")
        return SearchTreeType::Kdop;
    else
        return SearchTreeType::KdtreeInRadius;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
typename BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::CheckGap BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::ConvertCheckGap(const std::string& str)
{
    if(str == "NoCheck" || str == "no_check")
        return CheckGap::NoCheck;
    else if(str == "DirectCheck" || str == "direct_check")
        return CheckGap::DirectCheck;
    else if (str == "MappingCheck" || str == "mapping_check")
        return CheckGap::MappingCheck;
    else
        return CheckGap::MappingCheck;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Parameters BaseContactSearch<TDim, TNumNodes, TNumNodesMaster>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "allocation_size"                      : 1000,
        "bucket_size"                          : 4,
        "search_factor"                        : 3.5,
        "type_search"                          : "InRadius",
        "check_gap"                            : "MappingCheck",
        "condition_name"                       : "",
        "final_string"                         : "",
        "inverted_search"                      : false,
        "dynamic_search"                       : false,
        "static_check_movement"                : false,
        "predefined_master_slave"              : true,
        "id_name"                              : "",
        "consider_gap_threshold"               : false,
        "predict_correct_lagrange_multiplier"  : false,
        "debug_mode"                           : false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class BaseContactSearch<2, 2>;
template class BaseContactSearch<3, 3>;
template class BaseContactSearch<3, 4>;
template class BaseContactSearch<3, 3, 4>;
template class BaseContactSearch<3, 4, 3>;

}  // namespace Kratos.
