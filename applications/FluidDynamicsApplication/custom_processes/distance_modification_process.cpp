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
    const double FactorCoeff,
    const double DistanceThreshold,
    const bool CheckAtEachStep,
    const bool NegElemDeactivation,
    const bool RecoverOriginalDistance)
    : Process(), mrModelPart(rModelPart) {

    mFactorCoeff = FactorCoeff;
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
        "distance_threshold"                     : 0.01,
        "check_at_each_time_step"                : false,
        "deactivate_full_negative_elements"      : true,
        "recover_original_distance_at_each_step" : false
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mFactorCoeff = rParameters["distance_factor"].GetDouble();
    mDistanceThreshold = rParameters["distance_threshold"].GetDouble();
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
}

void DistanceModificationProcess::ExecuteInitialize() {

    KRATOS_TRY;

    // Required variables check
    const auto& r_node = *mrModelPart.NodesBegin();
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_H, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

    // Obtain NODAL_H values
    FindNodalHProcess NodalHCalculator(mrModelPart);
    NodalHCalculator.Execute();

    KRATOS_CATCH("");
}

void DistanceModificationProcess::ExecuteBeforeSolutionLoop() {

    KRATOS_TRY;

    unsigned int bad_cuts = 1;

    // Modify the nodal distance values until there is no bad intersections
    while (bad_cuts > 0) {
        bad_cuts = this->ModifyDistance();
        mDistanceThreshold /= mFactorCoeff;
    }

    // If proceeds (depending on the formulation), perform the deactivation
    // Deactivates the full negative elements and sets values in the full negative elements
    if (mNegElemDeactivation) {
        this->DeactivateFullNegativeElements();
    }

    KRATOS_CATCH("");
}

void DistanceModificationProcess::ExecuteInitializeSolutionStep() {

    if(mCheckAtEachStep == true) {
        DistanceModificationProcess::ExecuteBeforeSolutionLoop();
    }
}

void DistanceModificationProcess::ExecuteFinalizeSolutionStep() {

    if(mRecoverOriginalDistance == true) {
        RecoverOriginalDistance();
    }
}

/* Protected functions ****************************************************/

unsigned int DistanceModificationProcess::ModifyDistance() {

    ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

    // Distance modification
    // Case in where the original distance does not need to be preserved (e.g. CFD)
    if (mRecoverOriginalDistance == false) {
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rNodes.size()); ++k) {
            ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
            const double h = itNode->FastGetSolutionStepValue(NODAL_H);
            const double tol_d = mDistanceThreshold*h;
            double& d = itNode->FastGetSolutionStepValue(DISTANCE);

            // Modify the distance to avoid almost empty fluid elements
            if((d >= 0.0) && (d < tol_d))
                d = -0.001*tol_d;
        }
    }
    // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
    else {
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector<std::vector<unsigned int>> AuxModifiedDistancesIDs(NumThreads);
        std::vector<std::vector<double>> AuxModifiedDistancesValues(NumThreads);

        #pragma omp parallel shared(AuxModifiedDistancesIDs, AuxModifiedDistancesValues)
        {
            const int ThreadId = OpenMPUtils::ThisThread();             // Get the thread id
            std::vector<unsigned int>   LocalModifiedDistancesIDs;      // Local modified distances nodes id vector
            std::vector<double>      LocalModifiedDistancesValues;      // Local modified distances original values vector

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rNodes.size()); ++k) {
                ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
                const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                const double tol_d = mDistanceThreshold*h;
                double& d = itNode->FastGetSolutionStepValue(DISTANCE);

                if((d >= 0.0) && (d < tol_d)) {

                    // Store the original distance to be recovered at the end of the step
                    LocalModifiedDistancesIDs.push_back(d);
                    LocalModifiedDistancesValues.push_back(itNode->Id());

                    // Modify the distance to avoid almost empty fluid elements
                    d = -0.001*tol_d;
                }
            }

            AuxModifiedDistancesIDs[ThreadId] = LocalModifiedDistancesIDs;
            AuxModifiedDistancesValues[ThreadId] = LocalModifiedDistancesValues;
        }

        mModifiedDistancesIDs = AuxModifiedDistancesIDs;
        mModifiedDistancesValues = AuxModifiedDistancesValues;
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

    // Check if there still exist bad cuts
    int num_bad_cuts = 0;
    /* Note: I'm defining a temporary variable because 'num_bad_cuts'
    *  instead of writing directly into input argument 'bad_cuts'
    *  because using a reference variable in a reduction pragma does
    *  not compile in MSVC 2015 nor in clang-3.8 (Is it even allowed by omp?)
    */
    #pragma omp parallel for reduction(+ : num_bad_cuts)
    for (int k = 0; k < static_cast<int>(rElements.size()); ++k) {
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;

        ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
        GeometryType& rGeometry = itElement->GetGeometry();

        for (unsigned int iNode=0; iNode<rGeometry.size(); iNode++) {
            const double d = rGeometry[iNode].FastGetSolutionStepValue(DISTANCE);
            (d > 0.0) ? n_pos++ : n_neg++;
        }

        // The element is cut
        if((n_neg > 0) && (n_pos > 0)) {
            for(unsigned int iNode=0; iNode<rGeometry.size(); iNode++) {
                const Node<3> &rConstNode = rGeometry[iNode];
                const double h = rConstNode.GetValue(NODAL_H);
                const double tol_d = (mDistanceThreshold*mFactorCoeff)*h;
                const double d = rConstNode.FastGetSolutionStepValue(DISTANCE);

                if((d >= 0.0) && (d < tol_d)) {
                    num_bad_cuts++;
                    break;
                }
            }
        }
    }

    // Synchronize the number of bad cuts
    mrModelPart.GetCommunicator().SumAll(num_bad_cuts);

    return num_bad_cuts;
}

void DistanceModificationProcess::RecoverOriginalDistance() {
    #pragma omp parallel
    {
        const int ThreadId = OpenMPUtils::ThisThread();
        const std::vector<unsigned int> LocalModifiedDistancesIDs = mModifiedDistancesIDs[ThreadId];
        const std::vector<double> LocalModifiedDistancesValues = mModifiedDistancesValues[ThreadId];

        for(unsigned int i=0; i<LocalModifiedDistancesIDs.size(); ++i) {
            const unsigned int nodeId = LocalModifiedDistancesIDs[i];
            mrModelPart.GetNode(nodeId).FastGetSolutionStepValue(DISTANCE) = LocalModifiedDistancesValues[i];
        }
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

    // Empty the modified distance vectors
    mModifiedDistancesIDs.resize(0);
    mModifiedDistancesValues.resize(0);
    mModifiedDistancesIDs.shrink_to_fit();
    mModifiedDistancesValues.shrink_to_fit();
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
        GeometryType& rGeometry = itElement->GetGeometry();

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
