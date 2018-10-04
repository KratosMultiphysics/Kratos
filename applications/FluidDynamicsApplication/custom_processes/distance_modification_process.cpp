//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
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
    : Process(), mrModelPart(rModelPart) {

    mDistanceThreshold = DistanceThreshold;
    mCheckAtEachStep = CheckAtEachStep;
    mNegElemDeactivation = NegElemDeactivation;
    mRecoverOriginalDistance = RecoverOriginalDistance;
}

DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart) {

    Parameters default_parameters( R"(
    {
        "model_part_name"                        : "default_model_part_name",
        "distance_factor"                        : 2.0,
        "distance_threshold"                     : 0.001,
        "continuous_distance"                    : true,
        "check_at_each_time_step"                : false,
        "avoid_almost_empty_elements"            : true,
        "deactivate_full_negative_elements"      : true,
        "recover_original_distance_at_each_step" : false
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mIsModified = false;
    mDistanceThreshold = rParameters["distance_threshold"].GetDouble();
    mContinuousDistance = rParameters["continuous_distance"].GetBool();
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mAvoidAlmostEmptyElements = rParameters["avoid_almost_empty_elements"].GetBool();
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
}

void DistanceModificationProcess::ExecuteInitialize() {

    KRATOS_TRY;

    // Continuous distance field required variables check
    if (mContinuousDistance){
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_H, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
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
            // Compute NODAL_H (used for computing the distance tolerance)
            FindNodalHProcess<true> nodal_h_calculator(mrModelPart);
            nodal_h_calculator.Execute();
            // Modify the continuous distance field
            this->ModifyDistance();
        } else {
            // Modify the discontinuous distance field
            this->ModifyDiscontinuousDistance();
        }

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
            const double h = it_node->FastGetSolutionStepValue(NODAL_H);
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

            for (auto it_node = nodes_begin; it_node != nodes_end; ++it_node) {
                const double h = it_node->FastGetSolutionStepValue(NODAL_H);
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
                    r_elem_dist(i_node) = -tol_d;
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

            for (auto it_elem = elems_begin; it_elem != elems_end; ++it_elem){
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
                        r_elem_dist(i_node) = -tol_d;
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
}

void DistanceModificationProcess::RecoverDeactivationPreviousState(){
    // Activate again all the elements
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        it_elem->Set(ACTIVE,true);
    }

    // Free the negative DOFs that were fixed
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node){
        auto it_node = mrModelPart.NodesBegin() + i_node;
        if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
            // Fix the nodal DOFs
            it_node->Free(PRESSURE);
            it_node->Free(VELOCITY_X);
            it_node->Free(VELOCITY_Y);
            it_node->Free(VELOCITY_Z);
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
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        ModelPart::NodesContainerType::iterator it_node = rNodes.begin() + i_node;
        if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
            // Fix the nodal DOFs
            it_node->Fix(PRESSURE);
            it_node->Fix(VELOCITY_X);
            it_node->Fix(VELOCITY_Y);
            it_node->Fix(VELOCITY_Z);
            // Set to zero the nodal DOFs
            it_node->FastGetSolutionStepValue(PRESSURE) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
        }
    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.
