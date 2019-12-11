// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/mpc_contact_search_process.h"
#include "custom_master_slave_constraints/contact_master_slave_constraint.h"

/* Custom utilities */
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::MPCContactSearchProcess(
    ModelPart & rMainModelPart,
    Parameters ThisParameters
    ) : BaseType(rMainModelPart, ThisParameters)
{
    // If we are going to consider multple searchs
    const std::string& id_name = BaseType::mThisParameters["id_name"].GetString();

    // We get the contact model part
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = BaseType::mOptions.IsNot(BaseType::MULTIPLE_SEARCHS) ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub" + id_name);

    // Iterate in the constraints
    auto& r_constraints_array = r_sub_contact_model_part.MasterSlaveConstraints();
    const auto it_const_begin = r_constraints_array.begin();
    const int num_constraints = static_cast<int>(r_constraints_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_constraints; ++i)
        (it_const_begin + i)->Set(ACTIVE, false);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CheckContactModelParts()
{
    // Calling base class
    BaseType::CheckContactModelParts();

    // Iterate in the constraints
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = BaseType::mOptions.IsNot(BaseType::MULTIPLE_SEARCHS) ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+BaseType::mThisParameters["id_name"].GetString());
    auto& r_constraints_array = r_sub_contact_model_part.MasterSlaveConstraints();

    const SizeType total_number_constraints = BaseType::mrMainModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();

    std::vector<MasterSlaveConstraint::Pointer> auxiliar_constraints_vector;

    #pragma omp parallel
    {
        // Buffer for new constraints if created
        std::vector<MasterSlaveConstraint::Pointer> auxiliar_constraints_vector_buffer;

        #pragma omp for
        for(int i = 0; i < static_cast<int>(r_constraints_array.size()); ++i) {
            auto it_const = r_constraints_array.begin() + i;

            if (it_const->Is(MARKER)) {
                // Setting the flag to remove
                it_const->Set(TO_ERASE, true);

                // Creating new condition
                MasterSlaveConstraint::Pointer p_new_const = it_const->Clone(total_number_constraints + it_const->Id());
                auxiliar_constraints_vector_buffer.push_back(p_new_const);

                p_new_const->SetData(it_const->GetData()); // TODO: Remove when fixed on the core
                p_new_const->Set(Flags(*it_const));
                p_new_const->Set(MARKER, true);
            } else {
                // Setting the flag to mark
                it_const->Set(MARKER, true);
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(auxiliar_constraints_vector_buffer.begin(),auxiliar_constraints_vector_buffer.end(),back_inserter(auxiliar_constraints_vector));
        }
    }

    // Finally we add the new constraints to the model part
    r_sub_contact_model_part.RemoveMasterSlaveConstraints(TO_ERASE);
    // Reorder ids (in order to keep the ids consistent)
    for (int i = 0; i < static_cast<int>(auxiliar_constraints_vector.size()); ++i) {
        auxiliar_constraints_vector[i]->SetId(total_number_constraints + i + 1);
    }
    ModelPart::MasterSlaveConstraintContainerType aux_conds;
    aux_conds.GetContainer() = auxiliar_constraints_vector;
    r_sub_contact_model_part.AddMasterSlaveConstraints(aux_conds.begin(), aux_conds.end());

    // Unsetting TO_ERASE
    VariableUtils().SetFlag(TO_ERASE, false, r_contact_model_part.MasterSlaveConstraints());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
Condition::Pointer MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::AddPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId,
    GeometricalObject::Pointer pObjectSlave,
    const array_1d<double, 3>& rSlaveNormal,
    GeometricalObject::Pointer pObjectMaster,
    const array_1d<double, 3>& rMasterNormal,
    IndexMap::Pointer pIndexesPairs,
    Properties::Pointer pProperties
    )
{
    Condition::Pointer p_cond = BaseType::AddPairing(rComputingModelPart, rConditionId, pObjectSlave, rSlaveNormal, pObjectMaster, rMasterNormal, pIndexesPairs, pProperties);

    const bool is_frictional = BaseType::mrMainModelPart.Is(SLIP);
    const bool is_rigid = is_frictional ? false : BaseType::mrMainModelPart.Is(RIGID);

    // Creating constraint
    if (p_cond.get() != nullptr) {
        MasterSlaveConstraint::Pointer p_new_const = Kratos::make_shared<ContactMasterSlaveConstraint>(GetMaximumConstraintsIds() + 1);
        p_new_const->Set(ACTIVE);
        p_new_const->Initialize(rComputingModelPart.GetProcessInfo());
        rComputingModelPart.AddMasterSlaveConstraint(p_new_const);
        p_cond->SetValue(CONSTRAINT_POINTER, p_new_const);
        if (is_frictional) p_cond->Set(SLIP);
        if (is_rigid) p_cond->Set(RIGID);
        p_cond->Initialize();
    }

    return p_cond;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CleanModelPart(ModelPart& rModelPart)
{
    // Calling base class
    BaseType::CleanModelPart(rModelPart);

    // We clean the constraints
    auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
    VariableUtils().SetFlag(TO_ERASE, true, r_constraints_array);

    BaseType::mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline IndexType MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::GetMaximumConstraintsIds()
{
    auto& r_constraints_array = BaseType::mrMainModelPart.MasterSlaveConstraints();

    IndexType constraint_id = 0;
    for(IndexType i = 0; i < r_constraints_array.size(); ++i)  {
        auto it_const = r_constraints_array.begin() + i;
        const IndexType id = it_const->GetId();
        if (id > constraint_id)
            constraint_id = id;
    }

    return constraint_id;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::ResetContactOperators()
{
    // Calling the base class
    BaseType::ResetContactOperators();

    // We iterate over the master nodes
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = BaseType::mOptions.IsNot(BaseType::MULTIPLE_SEARCHS) ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub"+BaseType::mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();

    if (BaseType::mrMainModelPart.Is(MODIFIED)) { // It has been remeshed. We remove everything
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = r_nodes_array.begin() + i;
            if (it_node->Is(MASTER)) {
                IndexMap::Pointer p_indexes_pairs = it_node->GetValue(INDEX_MAP);

                if (p_indexes_pairs != nullptr) {
                    p_indexes_pairs->clear();
//                     p_indexes_pairs->reserve(mAllocationSize);
                }
            }
        }

        // We remove all the computing constraints
        const std::string sub_computing_model_part_name = "ComputingContactSub" + BaseType::mThisParameters["id_name"].GetString();
        ModelPart& r_computing_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("ComputingContact");
        ModelPart& r_sub_computing_contact_model_part = BaseType::mOptions.IsNot(BaseType::MULTIPLE_SEARCHS) ? r_computing_contact_model_part : r_computing_contact_model_part.GetSubModelPart(sub_computing_model_part_name);
        auto& r_computing_constraints_array = r_sub_computing_contact_model_part.MasterSlaveConstraints();
        const int num_computing_constraints = static_cast<int>(r_computing_constraints_array.size());

        #pragma omp parallel for
        for(int i = 0; i < num_computing_constraints; ++i) {
            auto it_const = r_computing_constraints_array.begin() + i;
            it_const->Set(TO_ERASE, true);
        }
    } else {
        // We iterate, but not in OMP
        for(IndexType i = 0; i < r_nodes_array.size(); ++i) {
            auto it_node = r_nodes_array.begin() + i;
            if (it_node->Is(MASTER)) {
                IndexMap::Pointer p_indexes_pairs = it_node->GetValue(INDEX_MAP);
                if (p_indexes_pairs != nullptr) {
                    // The vector with the ids to remove
                    std::vector<IndexType> inactive_constraints_ids;
                    for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                        MasterSlaveConstraint::Pointer p_const = BaseType::mrMainModelPart.pGetMasterSlaveConstraint(it_pair->second);
                        if (p_const->IsNot(ACTIVE)) {
                            p_const->Set(TO_ERASE, true);
                            inactive_constraints_ids.push_back(it_pair->first);
                        }
                    }
                    for (auto& i_to_remove : inactive_constraints_ids) {
                        p_indexes_pairs->RemoveId(inactive_constraints_ids[i_to_remove]);
                    }
                }
            }
        }
    }

    BaseType::mrMainModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
}


/***********************************************************************************/
/***********************************************************************************/

template class MPCContactSearchProcess<2, 2>;
template class MPCContactSearchProcess<3, 3>;
template class MPCContactSearchProcess<3, 4>;
template class MPCContactSearchProcess<3, 3, 4>;
template class MPCContactSearchProcess<3, 4, 3>;

}  // namespace Kratos.
