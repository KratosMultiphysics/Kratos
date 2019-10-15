//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "distance_modification_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    const double FactorCoeff, //TODO: Remove it (here for legacy reasons)
    const double DistanceThreshold,
    const bool CheckAtEachStep,
    const bool NegElemDeactivation,
    const bool RecoverOriginalDistance)
    : Process(),
      mrModelPart(rModelPart)
{
    // Member variables initialization
    mDistanceThreshold = DistanceThreshold;
    mCheckAtEachStep = CheckAtEachStep;
    mNegElemDeactivation = NegElemDeactivation;
    mRecoverOriginalDistance = RecoverOriginalDistance;

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    this->InitializeEmbeddedIsActive();
}

DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    this->InitializeEmbeddedIsActive();
}

DistanceModificationProcess::DistanceModificationProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    this->InitializeEmbeddedIsActive();
}

void DistanceModificationProcess::InitializeEmbeddedIsActive()
{
    for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node) {
        auto it_node = mrModelPart.NodesBegin() + i_node;
        it_node->SetValue(EMBEDDED_IS_ACTIVE, 0);
    }
}

void DistanceModificationProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    Parameters default_parameters( R"(
    {
        "model_part_name"                             : "",
        "distance_threshold"                          : 0.001,
        "continuous_distance"                         : true,
        "check_at_each_time_step"                     : true,
        "avoid_almost_empty_elements"                 : true,
        "deactivate_full_negative_elements"           : true,
        "recover_original_distance_at_each_step"      : false,
        "full_negative_elements_fixed_variables_list" : ["PRESSURE","VELOCITY"]
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mIsModified = false;
    mDistanceThreshold = rParameters["distance_threshold"].GetDouble();
    mContinuousDistance = rParameters["continuous_distance"].GetBool();
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mAvoidAlmostEmptyElements = rParameters["avoid_almost_empty_elements"].GetBool();
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
    if (mNegElemDeactivation) {
        this->CheckAndStoreVariablesList(rParameters["full_negative_elements_fixed_variables_list"].GetStringArray());
    }
}

void DistanceModificationProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void DistanceModificationProcess::ExecuteInitialize()
{
    KRATOS_TRY;

    if (mContinuousDistance){
        // Continuous distance modification historical variables check
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

        // Compute NODAL_H (used for computing the distance tolerance)
        FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable> nodal_h_calculator(mrModelPart);
        nodal_h_calculator.Execute();
    }

    KRATOS_CATCH("");
}

void DistanceModificationProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void DistanceModificationProcess::ExecuteInitializeSolutionStep() {

    if(!mIsModified){
        // Modify the nodal distance values to avoid bad intersections
        if (mContinuousDistance) {
            // Modify the continuous distance field
            this->ModifyDistance();
        } else {
            // Modify the discontinuous distance field
            this->ModifyDiscontinuousDistance();
        }
        mIsModified = true;

        // If proceeds (depending on the formulation), perform the deactivation
        // Deactivates the full negative elements and sets the inner values to 0
        if (mNegElemDeactivation) {
            this->DeactivateFullNegativeElements();
        }
    }
}

void DistanceModificationProcess::ExecuteFinalizeSolutionStep() {

    // Restore the initial state if the distance is checked at each time step
    if (mCheckAtEachStep){
        // Restore the is modified flag to false
        mIsModified = false;
        if (mNegElemDeactivation){
            // Recover the state previous to the element deactivation
            this->RecoverDeactivationPreviousState();
        }
    }

    // Recover the original distance values
    if(mRecoverOriginalDistance) {
        if (mContinuousDistance){
            this->RecoverOriginalDistance();
        } else {
            this->RecoverOriginalDiscontinuousDistance();
        }
    }
}

/* Protected functions ****************************************************/

void DistanceModificationProcess::ModifyDistance() {

    ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

    // Distance modification
    // Case in where the original distance does not need to be preserved (e.g. CFD)
    if (mRecoverOriginalDistance == false) {
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(r_nodes.size()); ++k) {
            auto it_node = r_nodes.begin() + k;
            const double h = it_node->GetValue(NODAL_H);
            const double tol_d = mDistanceThreshold*h;
            double& d = it_node->FastGetSolutionStepValue(DISTANCE);

            // Check if the distance values are close to zero
            // If proceeds, set the tolerance as distance value
            if(std::abs(d) < tol_d){
                if (d <= 0.0){
                    d = -tol_d;
                } else {
                    // If selected, avoid almost empty elements
                    if (mAvoidAlmostEmptyElements){
                        d = -tol_d;
                    } else {
                        d = tol_d;
                    }
                }
            }
        }
    }
    // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
    else {

        const int num_chunks = 2 * OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector partition_vec;
        OpenMPUtils::DivideInPartitions(r_nodes.size(),num_chunks,partition_vec);

        #pragma omp parallel for
        for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk)
        {
            auto nodes_begin = r_nodes.begin() + partition_vec[i_chunk];
            auto nodes_end = r_nodes.begin() + partition_vec[i_chunk + 1];

            // Auxiliar chunk arrays
            std::vector<unsigned int> aux_modified_distances_ids;
            std::vector<double> aux_modified_distance_values;

            for (auto it_node = nodes_begin; it_node < nodes_end; ++it_node) {
                const double h = it_node->GetValue(NODAL_H);
                const double tol_d = mDistanceThreshold*h;
                double &d = it_node->FastGetSolutionStepValue(DISTANCE);

                // Check if the distance values are close to zero
                // If proceeds, set the tolerance as distance value
                if(std::abs(d) < tol_d){

                    // Store the original distance to be recovered at the end of the step
                    aux_modified_distances_ids.push_back(it_node->Id());
                    aux_modified_distance_values.push_back(d);

                    if (d <= 0.0){
                        d = -tol_d;
                    } else {
                        // If selected, avoid almost empty elements
                        if (mAvoidAlmostEmptyElements){
                            d = -tol_d;
                        } else {
                            d = tol_d;
                        }
                    }
                }
            }

            // Save the auxiliar chunk arrays
            #pragma omp critical
            {
                mModifiedDistancesIDs.insert(mModifiedDistancesIDs.end(),aux_modified_distances_ids.begin(),aux_modified_distances_ids.end());
                mModifiedDistancesValues.insert(mModifiedDistancesValues.end(), aux_modified_distance_values.begin(), aux_modified_distance_values.end());
            }
        }
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

    // Update the TO_SPLIT flag
    this->SetContinuousDistanceToSplitFlag();
}

void DistanceModificationProcess::ModifyDiscontinuousDistance(){

    auto r_elems = mrModelPart.Elements();
    auto elems_begin = mrModelPart.ElementsBegin();
    const auto n_elems = mrModelPart.NumberOfElements();

    // Distance modification
    if (mRecoverOriginalDistance == false) {
        // Case in where the original distance does not need to be preserved (e.g. CFD)
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < static_cast<int>(n_elems); ++i_elem){
            auto it_elem = elems_begin + i_elem;

            // Compute the distance tolerance
            const double tol_d = mDistanceThreshold*(it_elem->GetGeometry()).Length();

            // Check if the distance values are close to zero
            Vector &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
            for (unsigned int i_node = 0; i_node < r_elem_dist.size(); ++i_node){
                if (std::abs(r_elem_dist(i_node)) < tol_d){
                    r_elem_dist(i_node) = r_elem_dist(i_node) > 0.0 ? tol_d : -tol_d;
                }
            }
        }
    } else {
        // Case in where the original distance needs to be kept to track the interface (e.g. FSI)

        const int num_chunks = 2 * OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector partition_vec;
        OpenMPUtils::DivideInPartitions(n_elems,num_chunks,partition_vec);

        #pragma omp parallel for
        for (int i_chunk = 0; i_chunk < num_chunks; ++i_chunk)
        {
            auto elems_begin = r_elems.begin() + partition_vec[i_chunk];
            auto elems_end = r_elems.begin() + partition_vec[i_chunk + 1];

            // Auxiliar chunk arrays
            std::vector<unsigned int> aux_modified_distances_ids;
            std::vector<Vector> aux_modified_elemental_distances;

            for (auto it_elem = elems_begin; it_elem < elems_end; ++it_elem){
                // Compute the distance tolerance
                const double tol_d = mDistanceThreshold * (it_elem->GetGeometry()).Length();

                bool is_saved = false;
                Vector &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
                for (unsigned int i_node = 0; i_node < r_elem_dist.size(); ++i_node){
                    if (std::abs(r_elem_dist(i_node)) < tol_d){
                        if (!is_saved){
                            aux_modified_distances_ids.push_back(it_elem->Id());
                            aux_modified_elemental_distances.push_back(r_elem_dist);
                        }
                        r_elem_dist(i_node) = r_elem_dist(i_node) > 0.0 ? tol_d : -tol_d;
                    }
                }
            }

            // Save the auxiliar chunk arrays
            #pragma omp critical
            {
                mModifiedDistancesIDs.insert(mModifiedDistancesIDs.end(),aux_modified_distances_ids.begin(),aux_modified_distances_ids.end());
                mModifiedElementalDistancesValues.insert(mModifiedElementalDistancesValues.end(),aux_modified_elemental_distances.begin(),aux_modified_elemental_distances.end());
            }
        }
    }

    // Update the TO_SPLIT flag
    this->SetDiscontinuousDistanceToSplitFlag();
}

void DistanceModificationProcess::RecoverDeactivationPreviousState(){
    // Activate again all the elements
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        it_elem->Set(ACTIVE,true);
    }
    if ((mDoubleVariablesList.size() > 0.0) || (mComponentVariablesList.size() > 0.0)){
        // Free the negative DOFs that were fixed
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node){
            auto it_node = mrModelPart.NodesBegin() + i_node;
            if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
                for (std::size_t i_var = 0; i_var < mDoubleVariablesList.size(); i_var++){
                    const auto& r_double_var = *mDoubleVariablesList[i_var];
                    // Free the nodal DOFs  that were fixed
                    it_node->Free(r_double_var);
                }
                for (std::size_t i_comp = 0; i_comp < mComponentVariablesList.size(); i_comp++){
                    const auto& r_component_var = *mComponentVariablesList[i_comp];
                    // Free the nodal DOFs that were fixed
                    it_node->Free(r_component_var);
                }
            }
        }
    }
}

void DistanceModificationProcess::RecoverOriginalDistance() {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mModifiedDistancesIDs.size()); ++i) {
        const auto node_id = mModifiedDistancesIDs[i];
        mrModelPart.GetNode(node_id).FastGetSolutionStepValue(DISTANCE) = mModifiedDistancesValues[i];
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

    // Empty the modified distance vectors
    mModifiedDistancesIDs.resize(0);
    mModifiedDistancesValues.resize(0);
    mModifiedDistancesIDs.shrink_to_fit();
    mModifiedDistancesValues.shrink_to_fit();

    // Restore the TO_SPLIT flag original status
    this->SetContinuousDistanceToSplitFlag();
}

void DistanceModificationProcess::RecoverOriginalDiscontinuousDistance() {
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mModifiedDistancesIDs.size()); ++i_elem) {
        const unsigned int elem_id = mModifiedDistancesIDs[i_elem];
        const auto elem_dist = mModifiedElementalDistancesValues[i_elem];
        mrModelPart.GetElement(elem_id).SetValue(ELEMENTAL_DISTANCES,elem_dist);
    }

    // Empty the modified distance vectors
    mModifiedDistancesIDs.resize(0);
    mModifiedElementalDistancesValues.resize(0);
    mModifiedDistancesIDs.shrink_to_fit();
    mModifiedElementalDistancesValues.shrink_to_fit();

    // Restore the TO_SPLIT flag original status
    this->SetDiscontinuousDistanceToSplitFlag();
}

void DistanceModificationProcess::DeactivateFullNegativeElements() {

    ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        ModelPart::NodesContainerType::iterator it_node = rNodes.begin() + i_node;
        it_node->SetValue(EMBEDDED_IS_ACTIVE, 0);
    }

    // Deactivate those elements whose negative distance nodes summation is equal to their number of nodes
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(rElements.size()); ++k){
        unsigned int n_neg = 0;
        ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
        auto& rGeometry = itElement->GetGeometry();

        // Check the distance function sign at the element nodes
        for (unsigned int i_node=0; i_node<rGeometry.size(); i_node++){
            if (rGeometry[i_node].FastGetSolutionStepValue(DISTANCE) < 0.0){
                n_neg++;
            }
        }

        (n_neg == rGeometry.size()) ? itElement->Set(ACTIVE, false) : itElement->Set(ACTIVE, true);

        // If the element is ACTIVE, all its nodes are active as well
        if (itElement->Is(ACTIVE)){
            for (unsigned int i_node = 0; i_node < rGeometry.size(); ++i_node){
                int& activation_index = rGeometry[i_node].GetValue(EMBEDDED_IS_ACTIVE);
                #pragma omp atomic
                activation_index += 1;
            }
        }
    }

    // Synchronize the EMBEDDED_IS_ACTIVE variable flag
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(EMBEDDED_IS_ACTIVE);

    // Set to zero and fix the DOFs in the remaining inactive nodes
    if ((mDoubleVariablesList.size() > 0.0) || (mComponentVariablesList.size() > 0.0)){
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
            auto it_node = rNodes.begin() + i_node;
            if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
                for (std::size_t i_var = 0; i_var < mDoubleVariablesList.size(); i_var++){
                    const auto& r_double_var = *mDoubleVariablesList[i_var];
                    // Fix the nodal DOFs
                    it_node->Fix(r_double_var);
                    // Set to zero the nodal DOFs
                    it_node->FastGetSolutionStepValue(r_double_var) = 0.0;
                }
                for (std::size_t i_comp = 0; i_comp < mComponentVariablesList.size(); i_comp++){
                    const auto& r_component_var = *mComponentVariablesList[i_comp];
                    // Fix the nodal DOFs
                    it_node->Fix(r_component_var);
                    // Set to zero the nodal DOFs
                    it_node->FastGetSolutionStepValue(r_component_var) = 0.0;
                }
            }
        }
    }
}

void DistanceModificationProcess::SetContinuousDistanceToSplitFlag()
{
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem) {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        auto &r_geom = it_elem->GetGeometry();
        std::vector<double> elem_dist;
        for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            elem_dist.push_back(r_geom[i_node].FastGetSolutionStepValue(DISTANCE));
        }
        this->SetElementToSplitFlag(*it_elem, elem_dist);
    }
}

void DistanceModificationProcess::SetDiscontinuousDistanceToSplitFlag()
{
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem) {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        const auto &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
        this->SetElementToSplitFlag(*it_elem, r_elem_dist);
    }
}

void DistanceModificationProcess::CheckAndStoreVariablesList(const std::vector<std::string>& rVariableStringArray)
{
    const auto& r_node = *mrModelPart.NodesBegin();
    for (std::size_t i_variable=0; i_variable < rVariableStringArray.size(); i_variable++){
        if (KratosComponents<Variable<double>>::Has(rVariableStringArray[i_variable])) {
            const auto& r_double_var  = KratosComponents<Variable<double>>::Get(rVariableStringArray[i_variable]);
            KRATOS_CHECK_DOF_IN_NODE(r_double_var, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_double_var, r_node)

            mDoubleVariablesList.push_back(&r_double_var);
        }
        else if (KratosComponents<ComponentType>::Has(rVariableStringArray[i_variable])){
            const auto& r_component_var  = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]);
            KRATOS_CHECK_DOF_IN_NODE(r_component_var, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_component_var, r_node)

            mComponentVariablesList.push_back(&r_component_var);
        }
        else if (KratosComponents<Variable<array_1d<double,3>>>::Has(rVariableStringArray[i_variable])){
            // Checking vector variable in nodal data
            const auto& r_vector_var = KratosComponents<Variable<array_1d<double,3>>>::Get(rVariableStringArray[i_variable]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_vector_var, r_node)

            // Checking and storing the component variables
            const auto& r_component_var_x = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_X");
            KRATOS_CHECK_DOF_IN_NODE(r_component_var_x, r_node);
            mComponentVariablesList.push_back(&r_component_var_x);

            const auto& r_component_var_y = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_Y");
            KRATOS_CHECK_DOF_IN_NODE(r_component_var_y, r_node);
            mComponentVariablesList.push_back(&r_component_var_y);

            if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
                const auto& r_component_var_z = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_Z");
                KRATOS_CHECK_DOF_IN_NODE(r_component_var_z, r_node);
                mComponentVariablesList.push_back(&r_component_var_z);
            }
        }
        else {
            KRATOS_ERROR << "The variable defined in the list is not a double variable nor a component variable. Given variable: " << rVariableStringArray[i_variable] << std::endl;
        }
    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.
