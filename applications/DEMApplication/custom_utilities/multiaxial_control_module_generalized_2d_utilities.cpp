//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Guillermo Casas
//

// Project includes
#include "custom_utilities/multiaxial_control_module_generalized_2d_utilities.hpp"

namespace Kratos
{

MultiaxialControlModuleGeneralized2DUtilities::MultiaxialControlModuleGeneralized2DUtilities(ModelPart& rDemModelPart,
                                ModelPart& rFemModelPart,
                                Parameters& rParameters):
                                mrDemModelPart(rDemModelPart),
                                mrFemModelPart(rFemModelPart) {
    KRATOS_TRY

    Parameters default_parameters( R"(
        {
            "dem_model_part_name"      : "SpheresPart",
            "fem_model_part_name"      : "RigidFacePart",
            "Parameters"    : {
                "control_module_delta_time": 1.0e-5,
                "perturbation_tolerance": 1.0e-2,
                "perturbation_period": 10,
                "max_reaction_rate_factor": 2.0,
                "stiffness_averaging_time_interval": 1.0,
                "velocity_averaging_time_interval": 1.0,
                "reaction_averaging_time_interval": 1.0,
                "output_interval": 1
            },
            "list_of_actuators" : []
        }  )" );

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mCMDeltaTime = rParameters["Parameters"]["control_module_delta_time"].GetDouble();
    mCMStep = 0;
    mStep = 0;
    mCMTime = 0.0;
    mKDeltaTime = rParameters["Parameters"]["stiffness_averaging_time_interval"].GetDouble();
    mKStep = 0;
    mKTime = 0.0;
    mActuatorCounter = 0;
    mPerturbationTolerance = rParameters["Parameters"]["perturbation_tolerance"].GetDouble();
    mPerturbationPeriod = rParameters["Parameters"]["perturbation_period"].GetInt();
    mMaxReactionCorrectionFraction = rParameters["Parameters"]["max_reaction_rate_factor"].GetDouble();
    const double stiffness_averaging_time_interval = rParameters["Parameters"]["stiffness_averaging_time_interval"].GetDouble();
    KRATOS_ERROR_IF(mCMDeltaTime > stiffness_averaging_time_interval) << "Control module delta time is bigger than stiffness_averaging_time_interval!"<< std::endl;
    mStiffnessAlpha = 1.0 - mCMDeltaTime / stiffness_averaging_time_interval; // NOTE: Regardless of using Ktime to calculate K, we still need to filter K.
    const double velocity_averaging_time_interval = rParameters["Parameters"]["velocity_averaging_time_interval"].GetDouble();
    KRATOS_ERROR_IF(mCMDeltaTime > velocity_averaging_time_interval) << "Control module delta time is bigger than velocity_averaging_time_interval!"<< std::endl;
    mVelocityAlpha = 1.0 - mCMDeltaTime / velocity_averaging_time_interval;
    const double reaction_averaging_time_interval = rParameters["Parameters"]["reaction_averaging_time_interval"].GetDouble();
    const double & dem_delta_time = mrDemModelPart.GetProcessInfo()[DELTA_TIME];
    KRATOS_ERROR_IF(dem_delta_time > reaction_averaging_time_interval) << "DEM delta time is bigger than reaction_averaging_time_interval!"<< std::endl;
    mReactionAlpha = 1.0 - dem_delta_time / reaction_averaging_time_interval;

    const unsigned int number_of_actuators = rParameters["list_of_actuators"].size();
    mVelocity.resize(number_of_actuators, false);
    mAcceleration.resize(number_of_actuators, false);
    mReactionStress.resize(number_of_actuators, false);
    mReactionStressOld.resize(number_of_actuators, false);
    mDisplacement.resize(number_of_actuators, false);
    mDisplacementOld.resize(number_of_actuators, false);
    mElasticReactionStress.resize(number_of_actuators, false);
    mStiffness.resize(number_of_actuators,number_of_actuators,false);
    noalias(mStiffness) = ZeroMatrix(number_of_actuators, number_of_actuators);
    mDeltaDisplacement.resize(number_of_actuators, number_of_actuators, false);
    noalias(mDeltaDisplacement) = ZeroMatrix(number_of_actuators, number_of_actuators);
    mDeltaReactionStress.resize(number_of_actuators, number_of_actuators, false);
    noalias(mDeltaReactionStress) = ZeroMatrix(number_of_actuators, number_of_actuators);
    ModelPart* p_sub_modelpart;
    for (unsigned int i = 0; i < number_of_actuators; i++) {
        Parameters this_actuator_parameters = rParameters["list_of_actuators"][i];
        const std::string& actuator_name = this_actuator_parameters["Parameters"]["actuator_name"].GetString();
        if(actuator_name == "Radial") {
            this_actuator_parameters.ValidateAndAssignDefaults(GetDefaultParametersForRadialActuator());
            const unsigned int number_of_fem_boundaries = this_actuator_parameters["list_of_fem_boundaries"].size();
            std::vector<ModelPart*> list_of_fem_submodelparts(number_of_fem_boundaries);
            std::vector<array_1d<double,3>> list_of_fem_outer_normals(number_of_fem_boundaries);
            for (unsigned int j = 0; j < number_of_fem_boundaries; j++) {
                p_sub_modelpart = &(mrFemModelPart.GetSubModelPart(this_actuator_parameters["list_of_fem_boundaries"][j]["model_part_name"].GetString()));
                list_of_fem_submodelparts[j] = p_sub_modelpart;
                AddTableToSubModelPart(i + 1, this_actuator_parameters["target_stress_table"], p_sub_modelpart);
            }
            mListsOfFEMSubModelPartsForEachActuator[actuator_name] = list_of_fem_submodelparts;
            mTargetStressTableIds[actuator_name] = i + 1;
        }
        else if (actuator_name == "RadialMultiDofs") {
            this_actuator_parameters.ValidateAndAssignDefaults(GetDefaultParametersForRadialMultiDofsActuator());
            const unsigned int number_of_fem_boundaries = this_actuator_parameters["list_of_fem_boundaries"].size();
            std::vector<ModelPart*> list_of_fem_submodelparts(number_of_fem_boundaries);
            std::vector<array_1d<double,3>> list_of_fem_outer_normals(number_of_fem_boundaries);
            for (unsigned int j = 0; j < number_of_fem_boundaries; j++) {
                p_sub_modelpart = &(mrFemModelPart.GetSubModelPart(this_actuator_parameters["list_of_fem_boundaries"][j]["model_part_name"].GetString()));
                list_of_fem_submodelparts[j] = p_sub_modelpart;
            }
            mListsOfFEMSubModelPartsForEachActuator[actuator_name] = list_of_fem_submodelparts;
            mMaxNodalVelocityForMultiDofs = this_actuator_parameters["Parameters"]["max_nodal_velocity"].GetDouble();
            mCompressionLengthForMultiDofs = this_actuator_parameters["Parameters"]["compression_length"].GetDouble();
        }
        else if (actuator_name == "Z") {
            this_actuator_parameters.ValidateAndAssignDefaults(GetDefaultParametersForZActuator());
            const unsigned int number_of_dem_boundaries = this_actuator_parameters["list_of_dem_boundaries"].size();
            std::vector<ModelPart*> list_of_dem_submodelparts(number_of_dem_boundaries);
            for (unsigned int j = 0; j < number_of_dem_boundaries; j++) {
                p_sub_modelpart = &(mrDemModelPart.GetSubModelPart(this_actuator_parameters["list_of_dem_boundaries"][j]["model_part_name"].GetString()));
                list_of_dem_submodelparts[j] = p_sub_modelpart;
                AddTableToSubModelPart(i + 1, this_actuator_parameters["target_stress_table"], p_sub_modelpart);
            }
            mListsOfDEMSubModelPartsForEachActuator[actuator_name] = list_of_dem_submodelparts;
            mTargetStressTableIds[actuator_name] = i + 1;
        }
        else if (actuator_name == "X" || actuator_name == "Y") {
            this_actuator_parameters.ValidateAndAssignDefaults(GetDefaultParametersForXOrYActuator());
            const unsigned int number_of_fem_boundaries = this_actuator_parameters["list_of_fem_boundaries"].size();
            std::vector<ModelPart*> list_of_fem_submodelparts(number_of_fem_boundaries);
            std::vector<array_1d<double,3>> list_of_fem_outer_normals(number_of_fem_boundaries);
            for (unsigned int j = 0; j < number_of_fem_boundaries; j++) {
                p_sub_modelpart = &(mrFemModelPart.GetSubModelPart(this_actuator_parameters["list_of_fem_boundaries"][j]["model_part_name"].GetString()));

                list_of_fem_submodelparts[j] = p_sub_modelpart;
                list_of_fem_outer_normals[j][0] = this_actuator_parameters["list_of_fem_boundaries"][j]["outer_normal"][0].GetDouble();
                list_of_fem_outer_normals[j][1] = this_actuator_parameters["list_of_fem_boundaries"][j]["outer_normal"][1].GetDouble();
                list_of_fem_outer_normals[j][2] = this_actuator_parameters["list_of_fem_boundaries"][j]["outer_normal"][2].GetDouble();
                // Unit normal vector
                double inverse_norm = 1.0/norm_2(list_of_fem_outer_normals[j]);
                list_of_fem_outer_normals[j][0] *= inverse_norm;
                list_of_fem_outer_normals[j][1] *= inverse_norm;
                list_of_fem_outer_normals[j][2] *= inverse_norm;
                AddTableToSubModelPart(i + 1, this_actuator_parameters["target_stress_table"], p_sub_modelpart);
            }
            mListsOfFEMSubModelPartsForEachActuator[actuator_name] = list_of_fem_submodelparts;
            mFEMOuterNormals[actuator_name] = list_of_fem_outer_normals;
            mTargetStressTableIds[actuator_name] = i + 1;
        }
        else {
            KRATOS_ERROR << "Unknown actuator type. Please check parameter 'actuator_name' of actuator in position " << i <<"."<<std::endl;
        }

        const double compression_length = this_actuator_parameters["Parameters"]["compression_length"].GetDouble();
        const double initial_velocity = this_actuator_parameters["Parameters"]["initial_velocity"].GetDouble();
        const double stiffness = this_actuator_parameters["Parameters"]["young_modulus"].GetDouble()/compression_length; // mStiffness is actually a stiffness over an area
        mVelocity[i] = initial_velocity;
        mStiffness(i,i) = stiffness;
        mAcceleration[i] = 0.0;
        mReactionStress[i] = 0.0;
        mReactionStressOld[i] = 0.0;
        mDisplacement[i] = 0.0;
        mDisplacementOld[i] = 0.0;
        mElasticReactionStress[i] = 0.0;
        mVectorOfActuatorNames.push_back(actuator_name);
    }

    mNormOfInitiallyEstimatedStiffness = 0.0;
    for(unsigned int i = 0; i < mStiffness.size1(); i++) {
        mNormOfInitiallyEstimatedStiffness += mStiffness(i,i) * mStiffness(i,i);
    }
    mNormOfInitiallyEstimatedStiffness = std::sqrt(mNormOfInitiallyEstimatedStiffness);

    // Initialize Variables
    mrDemModelPart.GetProcessInfo().SetValue(TARGET_STRESS_Z,0.0);
    array_1d<double,3> zero_vector = ZeroVector(3);
    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        const std::vector<ModelPart*>& DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
        // Iterate through all FEMBoundaries
        for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
            ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int j = 0; j<number_of_nodes; j++) {
                ModelPart::NodesContainerType::iterator it = it_begin + j;
                it->SetValue(TARGET_STRESS,zero_vector);
                it->SetValue(REACTION_STRESS,zero_vector);
                it->SetValue(SMOOTHED_REACTION_STRESS,zero_vector);
                it->SetValue(ELASTIC_REACTION_STRESS,zero_vector);
                it->SetValue(SMOOTHED_ELASTIC_REACTION_STRESS,zero_vector);
                it->SetValue(LOADING_VELOCITY,zero_vector);
            }
        }
        // Iterate through all DEMBoundaries
        for (unsigned int i = 0; i < DEMSubModelPartList.size(); i++) {
            ModelPart& rSubModelPart = *(DEMSubModelPartList[i]);
            // Iterate through nodes of DEM boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int j = 0; j<number_of_nodes; j++) {
                ModelPart::NodesContainerType::iterator it = it_begin + j;
                it->SetValue(TARGET_STRESS,zero_vector);
                it->SetValue(REACTION_STRESS,zero_vector);
                it->SetValue(SMOOTHED_REACTION_STRESS,zero_vector);
                it->SetValue(ELASTIC_REACTION_STRESS,zero_vector);
                it->SetValue(SMOOTHED_ELASTIC_REACTION_STRESS,zero_vector);
                it->SetValue(LOADING_VELOCITY,zero_vector);
            }
        }
    }

    CalculateMaximumInputStressVariationRate();

    KRATOS_CATCH("");
}

// Before DEM solution
void MultiaxialControlModuleGeneralized2DUtilities::ExecuteInitialize() {
    KRATOS_TRY;

    // Iterate through all on plane actuators (no radial or Z actuators) to set velocities of all nodes to 0.0
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "X" || actuator_name == "Y") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_velocity) = ZeroVector(3);
                }
            }
        } else if (actuator_name == "RadialMultiDofs") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->SetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY, 0.0);
                }
            }
        }
    }

    SetProvidedInitialVelocityToTheControlledBoundaries();

    KRATOS_CATCH("");
}

// Before DEM solution
void MultiaxialControlModuleGeneralized2DUtilities::ExecuteInitializeSolutionStep() {
    KRATOS_TRY;
    const double current_time = mrDemModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrDemModelPart.GetProcessInfo()[DELTA_TIME];
    mStep++;

    // Update velocities
    if (current_time > (mCMTime + 0.5 * delta_time)) {

        // Advance CM time
        mCMTime += mCMDeltaTime;
        mCMStep += 1;

        const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

        Vector next_target_stress(number_of_actuators);
        noalias(next_target_stress) = ZeroVector(number_of_actuators);

        // Iterate through all actuators
        for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
            const std::string& actuator_name = mVectorOfActuatorNames[ind];
            if(actuator_name == "RadialMultiDofs") { continue; }
            const std::vector<ModelPart*>& FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
            const std::vector<ModelPart*>& DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
            unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
            if (actuator_name == "Z") {
                TableType::Pointer pDEMTargetStressTable = (*(DEMSubModelPartList[0])).pGetTable(target_stress_table_id);
                next_target_stress[ind] = pDEMTargetStressTable->GetValue(mCMTime);
            } else {
                TableType::Pointer pFEMTargetStressTable = (*(FEMSubModelPartList[0])).pGetTable(target_stress_table_id);
                next_target_stress[ind] = pFEMTargetStressTable->GetValue(mCMTime);
            }
        }

        Vector target_stress_perturbation(number_of_actuators);
        noalias(target_stress_perturbation) = GetPerturbations(next_target_stress,mCMTime);
        noalias(next_target_stress) += target_stress_perturbation;

        // Calculate velocity
        CalculateVelocity(next_target_stress, current_time);

        // TODO: possible CM enhancement
            // We could calculate accelerations and limit them instead of calculating velocities and limiting their variation
        // CalculateAcceleration(next_target_stress, current_time);
    }

    // Move Actuators

    // TODO: possible CM enhancement
    // CalculateVelocity(next_target_stress, current_time);

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                r_velocity[0] = mVelocity[ind] * cos_theta;
                r_velocity[1] = mVelocity[ind] * sin_theta;
                r_velocity[2] = 0.0;
                noalias(r_delta_displacement) = r_velocity * delta_time;
                noalias(r_displacement) += r_delta_displacement;
                noalias(it->Coordinates()) = it->GetInitialPosition().Coordinates() + r_displacement;
            }
        } else if (actuator_name == "RadialMultiDofs") {
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                const double& cm_computed_velocity = it->GetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY);
                r_velocity[0] = cm_computed_velocity * cos_theta;
                r_velocity[1] = cm_computed_velocity * sin_theta;
                r_velocity[2] = 0.0;
                noalias(r_delta_displacement) = r_velocity * delta_time;
                noalias(r_displacement) += r_delta_displacement;
                noalias(it->Coordinates()) = it->GetInitialPosition().Coordinates() + r_displacement;
            }
        } else if (actuator_name == "Z") {
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] += mVelocity[ind]*delta_time/1.0;
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_velocity) = mVelocity[ind] * mFEMOuterNormals[actuator_name][i];
                    r_velocity[2] = 0.0;
                    noalias(r_delta_displacement) = r_velocity * delta_time;
                    noalias(r_displacement) += r_delta_displacement;
                    noalias(it->Coordinates()) = it->GetInitialPosition().Coordinates() + r_displacement;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

// After DEM solution
void MultiaxialControlModuleGeneralized2DUtilities::ExecuteFinalizeSolutionStep() {
    KRATOS_TRY
    const double current_time = mrDemModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrDemModelPart.GetProcessInfo()[DELTA_TIME];
    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

    // Update ReactionStresses
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "RadialMultiDofs") {
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& current_reaction_stress = it->GetValue(REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_reaction_stress = it->GetValue(SMOOTHED_REACTION_STRESS);
                    array_1d<double,3>& current_elastic_stress = it->GetValue(ELASTIC_REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_elastic_reaction_stress = it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS);

                    array_1d<double,3>& r_contact_force = it->FastGetSolutionStepValue(CONTACT_FORCES);
                    array_1d<double,3>& r_elastic_contact_force = it->FastGetSolutionStepValue(ELASTIC_FORCES);
                    current_reaction_stress = r_contact_force / it->FastGetSolutionStepValue(DEM_NODAL_AREA);
                    r_smoothed_reaction_stress = (1.0 - mReactionAlpha) * current_reaction_stress + mReactionAlpha * r_smoothed_reaction_stress;
                    current_elastic_stress = r_elastic_contact_force / it->FastGetSolutionStepValue(DEM_NODAL_AREA);
                    r_smoothed_elastic_reaction_stress = (1.0 - mReactionAlpha) * current_elastic_stress + mReactionAlpha * r_smoothed_elastic_reaction_stress;
                }
            }
        }
    }

    Vector reaction_stress_estimated(number_of_actuators);
    noalias(reaction_stress_estimated) = MeasureReactionStress(CONTACT_FORCES); // Total forces
    noalias(mReactionStress) = (1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress;
    // noalias(mReactionStress) = 1.0/(1.0 - std::pow(mReactionAlpha,mStep)) *
    //                         ((1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress);
    Vector elastic_reaction_stress_estimated(number_of_actuators);
    noalias(elastic_reaction_stress_estimated) = MeasureReactionStress(ELASTIC_FORCES); // Total forces - viscous forces
    noalias(mElasticReactionStress) = (1.0 - mReactionAlpha) * elastic_reaction_stress_estimated + mReactionAlpha * mElasticReactionStress;

    noalias(mDisplacement) += mVelocity * mCMDeltaTime;

    // Update Stiffness matrix
    if (current_time > (mKTime - 0.5 * delta_time)) {
        // Update K if DeltaDisplacement is invertible
        CalculateStiffness();
    }

    // Print results

    // Iterate through all on plane actuators to set variables to 0.0
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "X" || actuator_name == "Y") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_target_stress = it->GetValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->GetValue(REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_reaction_stress = it->GetValue(SMOOTHED_REACTION_STRESS);
                    array_1d<double,3>& r_elastic_reaction_stress = it->GetValue(ELASTIC_REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_elastic_reaction_stress = it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->GetValue(LOADING_VELOCITY);
                    noalias(r_target_stress) = ZeroVector(3);
                    noalias(r_reaction_stress) = ZeroVector(3);
                    noalias(r_smoothed_reaction_stress) = ZeroVector(3);
                    noalias(r_elastic_reaction_stress) = ZeroVector(3);
                    noalias(r_smoothed_elastic_reaction_stress) = ZeroVector(3);
                    noalias(r_loading_velocity) = ZeroVector(3);
                }
            }
        }
    }

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        const std::vector<ModelPart*>& DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
        unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(FEMSubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
            double current_target_stress = TargetStressTable->GetValue(current_time);
            #pragma omp parallel for
            for(int i = 0; i<number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                it->GetValue(TARGET_STRESS_X) = current_target_stress * cos_theta;
                it->GetValue(TARGET_STRESS_Y) = current_target_stress * sin_theta;
                it->GetValue(REACTION_STRESS_X) = reaction_stress_estimated[ind] * cos_theta;
                it->GetValue(REACTION_STRESS_Y) = reaction_stress_estimated[ind] * sin_theta;
                it->GetValue(SMOOTHED_REACTION_STRESS_X) = mReactionStress[ind] * cos_theta;
                it->GetValue(SMOOTHED_REACTION_STRESS_Y) = mReactionStress[ind] * sin_theta;
                it->GetValue(ELASTIC_REACTION_STRESS_X) = elastic_reaction_stress_estimated[ind] * cos_theta;
                it->GetValue(ELASTIC_REACTION_STRESS_Y) = elastic_reaction_stress_estimated[ind] * sin_theta;
                it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS_X) = mElasticReactionStress[ind] * cos_theta;
                it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS_Y) = mElasticReactionStress[ind] * sin_theta;
                it->GetValue(LOADING_VELOCITY_X) = mVelocity[ind] * cos_theta;
                it->GetValue(LOADING_VELOCITY_Y) = mVelocity[ind] * sin_theta;
            }
        } else if (actuator_name == "RadialMultiDofs") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(FEMSubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                it->GetValue(LOADING_VELOCITY_X) = mVelocity[ind] * cos_theta;
                it->GetValue(LOADING_VELOCITY_Y) = mVelocity[ind] * sin_theta;
            }
        } else if (actuator_name == "Z") {
            // Iterate through all DEMBoundaries
            for (unsigned int i = 0; i < DEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(DEMSubModelPartList[i]);
                // Iterate through nodes of DEM boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
                double current_target_stress = TargetStressTable->GetValue(current_time);
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->GetValue(TARGET_STRESS_Z) = current_target_stress;
                    it->GetValue(REACTION_STRESS_Z) = reaction_stress_estimated[ind];
                    it->GetValue(SMOOTHED_REACTION_STRESS_Z) = mReactionStress[ind];
                    it->GetValue(ELASTIC_REACTION_STRESS_Z) = elastic_reaction_stress_estimated[ind];
                    it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS_Z) = mElasticReactionStress[ind];
                    it->GetValue(LOADING_VELOCITY_Z) = mVelocity[ind];
                }
                mrDemModelPart.GetProcessInfo()[TARGET_STRESS_Z] = std::abs(current_target_stress);
            }
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
                double current_target_stress = TargetStressTable->GetValue(current_time);
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_target_stress = it->GetValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->GetValue(REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_reaction_stress = it->GetValue(SMOOTHED_REACTION_STRESS);
                    array_1d<double,3>& r_elastic_reaction_stress = it->GetValue(ELASTIC_REACTION_STRESS);
                    array_1d<double,3>& r_smoothed_elastic_reaction_stress = it->GetValue(SMOOTHED_ELASTIC_REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->GetValue(LOADING_VELOCITY);
                    noalias(r_target_stress) += current_target_stress * mFEMOuterNormals[actuator_name][i];
                    noalias(r_reaction_stress) += reaction_stress_estimated[ind] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_smoothed_reaction_stress) += mReactionStress[ind] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_elastic_reaction_stress) += elastic_reaction_stress_estimated[ind] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_smoothed_elastic_reaction_stress) += mElasticReactionStress[ind] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_loading_velocity) += mVelocity[ind] * mFEMOuterNormals[actuator_name][i];
                }
            }
        }
    }
    KRATOS_CATCH("");
}

//***************************************************************************************************************

Vector MultiaxialControlModuleGeneralized2DUtilities::MeasureReactionStress(const Variable<array_1d<double,3>>& rVariable) {
    KRATOS_TRY
    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();
    Vector reaction_stress(number_of_actuators);
    noalias(reaction_stress) = ZeroVector(number_of_actuators);

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        const std::vector<ModelPart*>& FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        const std::vector<ModelPart*>& DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
        double face_area = 0.0;
        double face_reaction = 0.0;
        if (actuator_name == "Radial") {
            // Calculate face_area
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through conditions of FEM boundary
                const int number_of_conditions = static_cast<int>(rSubModelPart.Conditions().size());
                ModelPart::ConditionsContainerType::iterator con_begin = rSubModelPart.ConditionsBegin();
                #pragma omp parallel for reduction(+:face_area)
                for(int j = 0; j < number_of_conditions; j++) {
                    ModelPart::ConditionsContainerType::iterator itCond = con_begin + j;
                    face_area += itCond->GetGeometry().Area();
                }
            }
            // Calculate face_reaction
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for reduction(+:face_reaction)
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(rVariable);
                    // Unit normal vector pointing outwards
                    array_1d<double,3> radial_normal;
                    radial_normal[0] = it->X();
                    radial_normal[1] = it->Y();
                    radial_normal[2] = 0.0;
                    double inv_norm = 1.0/norm_2(radial_normal);
                    radial_normal[0] *= inv_norm;
                    radial_normal[1] *= inv_norm;
                    face_reaction += -inner_prod(r_force,radial_normal);
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[ind] = face_reaction/face_area;
            } else {
                reaction_stress[ind] = 0.0;
            }
        } else if (actuator_name == "RadialMultiDofs") {


        } else if (actuator_name == "Z") {
            // Calculate face_area
            // Iterate through all DEMBoundaries
            for (unsigned int i = 0; i < DEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(DEMSubModelPartList[i]);
                // Iterate through elements of DEM boundary
                ModelPart::ElementsContainerType& rElements = rSubModelPart.GetCommunicator().LocalMesh().Elements();
                #pragma omp parallel for reduction(+:face_area)
                for (int j = 0; j < (int)rElements.size(); j++) {
                    ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + j;
                    Element* p_element = ptr_itElem->get();
                    SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
                    const double radius = pDemElem->GetRadius();
                    face_area += Globals::Pi*radius*radius;
                }
            }
            // Calculate face_reaction
            // Iterate through all DEMBoundaries
            for (unsigned int i = 0; i < DEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(DEMSubModelPartList[i]);
                // Iterate through elements of DEM boundary
                ModelPart::ElementsContainerType& rElements = rSubModelPart.GetCommunicator().LocalMesh().Elements();
                #pragma omp parallel for reduction(+:face_reaction)
                for (int j = 0; j < (int)rElements.size(); j++) {
                    ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + j;
                    Element* p_element = ptr_itElem->get();
                    SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
                    BoundedMatrix<double, 3, 3> stress_tensor = ZeroMatrix(3,3);
                    noalias(stress_tensor) = (*(pDemElem->mSymmStressTensor));
                    const double radius = pDemElem->GetRadius();
                    face_reaction += stress_tensor(2,2) * Globals::Pi*radius*radius;
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[ind] = face_reaction/face_area;
            } else {
                reaction_stress[ind] = 0.0;
            }
        } else {
            // Calculate face_area
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through conditions of FEM boundary
                const int number_of_conditions = static_cast<int>(rSubModelPart.Conditions().size());
                ModelPart::ConditionsContainerType::iterator con_begin = rSubModelPart.ConditionsBegin();
                #pragma omp parallel for reduction(+:face_area)
                for(int j = 0; j < number_of_conditions; j++) {
                    ModelPart::ConditionsContainerType::iterator itCond = con_begin + j;
                    face_area += itCond->GetGeometry().Area();
                }
            }
            // Calculate face_reaction
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for reduction(+:face_reaction)
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(rVariable);
                    face_reaction += -inner_prod(r_force,mFEMOuterNormals[actuator_name][i]);
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[ind] = face_reaction/face_area;
            } else {
                reaction_stress[ind] = 0.0;
            }
        }
    }

    return reaction_stress;
    KRATOS_CATCH("");
}

//***************************************************************************************************************

Vector MultiaxialControlModuleGeneralized2DUtilities::GetPerturbations(const Vector& rTargetStress, const double& rTime) {
    KRATOS_TRY
    const unsigned int number_of_actuators = rTargetStress.size();
    Vector stress_perturbation(number_of_actuators);
    noalias(stress_perturbation) = ZeroVector(number_of_actuators);

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        if (actuator_name == "Z") {
            stress_perturbation[ind] = 0.0;
        } else {
            double amplitude = rTargetStress[ind] * mPerturbationTolerance;
            double omega = 2.0 * Globals::Pi / (mPerturbationPeriod * mCMDeltaTime);
            double phi = ind * 2.0 * Globals::Pi / number_of_actuators;
            stress_perturbation[ind] = amplitude * std::sin(omega * rTime + phi);
        }
    }
    return stress_perturbation;
    KRATOS_CATCH("");
}

//***************************************************************************************************************

double MultiaxialControlModuleGeneralized2DUtilities::GetConditionNumber(const Matrix& rInputMatrix, const Matrix& rInvertedMatrix) {
    KRATOS_TRY
    // Find the condition number to define is inverse is OK
    const double input_matrix_norm = norm_frobenius(rInputMatrix);
    const double inverted_matrix_norm = norm_frobenius(rInvertedMatrix);

    const double cond_number = input_matrix_norm * inverted_matrix_norm ;

    return cond_number;
    KRATOS_CATCH("");
}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::CalculateVelocity(const Vector& r_next_target_stress, const double& r_current_time) {
    KRATOS_TRY
    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

    Vector delta_target_stress(number_of_actuators);
    noalias(delta_target_stress) = r_next_target_stress-mReactionStress;
    Vector velocity_estimated(number_of_actuators);

    Matrix k_inverse(number_of_actuators, number_of_actuators);
    double k_det = 0.0;
    MathUtils<double>::InvertMatrix(mStiffness, k_inverse, k_det, -1.0);
    const bool is_k_invertible = MathUtils<double>::CheckConditionNumber(mStiffness,
                                                    k_inverse, std::numeric_limits<double>::epsilon(),
                                                    false);
    const double k_condition_number = GetConditionNumber(mStiffness,k_inverse);

    Vector velocity_perturbation(number_of_actuators);
    noalias(velocity_perturbation) = GetPerturbations(mVelocity,r_current_time);
    if (is_k_invertible == false || std::isnan(k_condition_number)) {
        noalias(velocity_estimated) = mVelocity + velocity_perturbation;
        //KRATOS_WARNING("DEM") << "Stiffness matrix is not invertible. Keeping loading velocity constant" << std::endl;
    } else {
        noalias(velocity_estimated) = prod(k_inverse,delta_target_stress)/mCMDeltaTime;
    }

    double norm_stiffness = 0.0;
    for(unsigned int i = 0; i < mStiffness.size1(); i++) {
        norm_stiffness += mStiffness(i,i)*mStiffness(i,i);
    }
    norm_stiffness = std::sqrt(norm_stiffness);
    const double max_allowed_velocity = mMaxReactionCorrectionFraction * mMaximumInputStressVariationRate / mNormOfInitiallyEstimatedStiffness;
    const double norm_velocity = norm_2(velocity_estimated);
    if (norm_velocity > max_allowed_velocity) { //TODO: IP, why are you comparing with the norm and not each component separately?
        for (unsigned int i = 0; i < velocity_estimated.size(); i++) {
            velocity_estimated[i] = max_allowed_velocity/norm_velocity * velocity_estimated[i];
        }
    }
    noalias(mVelocity) = (1.0 - mVelocityAlpha) * velocity_estimated + mVelocityAlpha * mVelocity;

    // TODO: possible CM enhancement
    // noalias(mVelocity) += mAcceleration * mrDemModelPart.GetProcessInfo()[DELTA_TIME];
    // const double norm_velocity = norm_2(mVelocity);
    // if (norm_velocity > max_allowed_velocity) {
    //     for (unsigned int i = 0; i < mVelocity.size(); i++) {
    //         mVelocity[i] = max_allowed_velocity/norm_velocity * mVelocity[i];
    //     }
    // }

    for (unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        if (actuator_name == "RadialMultiDofs") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for (int i = 0; i < number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double& r_nodal_next_target_stress = it->GetValue(RADIAL_NORMAL_STRESS_COMPONENT); // This is set by an external tool (can be a value projected from another modelpart)
                const array_1d<double, 3>& r_nodal_reaction_stress_vector = it->GetValue(SMOOTHED_REACTION_STRESS);
                // Unit normal vector pointing outwards
                array_1d<double,3> radial_normal;
                radial_normal[0] = it->X();
                radial_normal[1] = it->Y();
                radial_normal[2] = 0.0;
                double inv_norm = 1.0/norm_2(radial_normal);
                radial_normal[0] *= inv_norm;
                radial_normal[1] *= inv_norm;
                const double nodal_reaction_stress = r_nodal_reaction_stress_vector[0] * radial_normal[0] + r_nodal_reaction_stress_vector[1] * radial_normal[1];

                const double nodal_delta_target_stress = r_nodal_next_target_stress + nodal_reaction_stress;

                if (std::abs(nodal_reaction_stress) < std::numeric_limits<double>::epsilon()) {
                    double nodal_smoothed_nodal_velocity = (1.0 - mVelocityAlpha) * -1.0 * mMaxNodalVelocityForMultiDofs + mVelocityAlpha * it->GetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY);
                    it->SetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY, nodal_smoothed_nodal_velocity);
                }
                else { //if (nodal_delta_target_stress < 0.0 && nodal_reaction_stress > 0.0) {
                    double nodal_velocity_estimated = nodal_delta_target_stress * mCompressionLengthForMultiDofs / mNormOfInitiallyEstimatedStiffness / mCMDeltaTime; //TODO: norm_stiffness?
                    if (std::abs(nodal_velocity_estimated) > mMaxNodalVelocityForMultiDofs) {
                        nodal_velocity_estimated *= mMaxNodalVelocityForMultiDofs / std::abs(nodal_velocity_estimated);
                    }
                    double nodal_smoothed_nodal_velocity = (1.0 - mVelocityAlpha) * nodal_velocity_estimated + mVelocityAlpha * it->GetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY);
                    it->SetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY, nodal_smoothed_nodal_velocity);
                }
                /*else {
                    it->SetValue(SMOOTHED_SCALAR_RADIAL_VELOCITY, 0.0);
                }*/
            }
        }
    }
    KRATOS_CATCH("");
}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::CalculateAcceleration(const Vector& r_next_target_stress, const double& r_current_time) {
    KRATOS_TRY
    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

    Vector delta_target_stress(number_of_actuators);
    noalias(delta_target_stress) = r_next_target_stress-mReactionStress;

    Matrix k_inverse(number_of_actuators,number_of_actuators);
    double k_det = 0.0;
    MathUtils<double>::InvertMatrix(mStiffness, k_inverse, k_det, -1.0);
    const bool is_k_invertible = MathUtils<double>::CheckConditionNumber(mStiffness,
                                                    k_inverse, std::numeric_limits<double>::epsilon(),
                                                    false);
    const double k_condition_number = GetConditionNumber(mStiffness,k_inverse);

    Vector acceleration_perturbation(number_of_actuators);
    noalias(acceleration_perturbation) = GetPerturbations(mAcceleration,r_current_time);
    if (is_k_invertible == false || std::isnan(k_condition_number)) {
        noalias(mAcceleration) = mAcceleration + acceleration_perturbation;
        std::cout << "Stiffness matrix is not invertible. Keeping acceleration constant" << std::endl;
    } else {
        noalias(mAcceleration) = 2.0 / (mCMDeltaTime * mCMDeltaTime) *
                                 (prod(k_inverse, delta_target_stress) - mVelocity * mCMDeltaTime);
    }

    // TODO: possible CM enhancement
    // Investigate more on how to limit acceleration
    // Limit acceleration
    double norm_stiffness = 0.0;
    for(unsigned int i = 0; i < mStiffness.size1(); i++) {
        norm_stiffness += mStiffness(i,i)*mStiffness(i,i);
    }
    norm_stiffness = std::sqrt(norm_stiffness);
    const double max_allowed_acceleration = mMaxReactionCorrectionFraction * mMaximumInputStressVariationRate / (norm_stiffness * mCMDeltaTime);
    const double norm_acceleration = norm_2(mAcceleration);
    if (norm_acceleration > max_allowed_acceleration) {
        for (unsigned int i = 0; i < mAcceleration.size(); i++) {
            mAcceleration[i] = max_allowed_acceleration/norm_acceleration * mAcceleration[i];
        }
    }
    KRATOS_CATCH("");
}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::CalculateStiffness() {
    KRATOS_TRY
    mKTime += mKDeltaTime;
    mKStep++;

    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

    Vector delta_reaction_stress(number_of_actuators);
    noalias(delta_reaction_stress) = mReactionStress - mReactionStressOld;
    noalias(mReactionStressOld) = mReactionStress;

    Vector delta_displacement(number_of_actuators);
    noalias(delta_displacement) = mDisplacement - mDisplacementOld;
    noalias(mDisplacementOld) = mDisplacement;

    for (unsigned int i = 0; i < number_of_actuators; i++) {
        mDeltaDisplacement(i,mActuatorCounter) = delta_displacement[i];
        mDeltaReactionStress(i,mActuatorCounter) = delta_reaction_stress[i];
    }

    Matrix k_estimated(number_of_actuators,number_of_actuators);

    if (mKStep > number_of_actuators-1) {
        Matrix delta_displacement_inverse(number_of_actuators, number_of_actuators);
        double delta_displacement_det = 0.0;
        MathUtils<double>::InvertMatrix(mDeltaDisplacement, delta_displacement_inverse, delta_displacement_det,-1.0);
        const bool is_delta_displacement_invertible = MathUtils<double>::CheckConditionNumber(mDeltaDisplacement,
                                                        delta_displacement_inverse, 1.0e-10,
                                                        false);
        const double delta_displacement_condition_number = GetConditionNumber(mDeltaDisplacement,delta_displacement_inverse);

        if (is_delta_displacement_invertible == false || std::isnan(delta_displacement_condition_number) ) {
            noalias(k_estimated) = mStiffness;
            //KRATOS_WARNING("DEM") << "Delta displacement matrix is not invertible. Keeping stiffness matrix constant" << std::endl;
        } else {
            noalias(k_estimated) = prod(mDeltaReactionStress,delta_displacement_inverse);
        }
        noalias(mStiffness) = (1.0 - mStiffnessAlpha) * k_estimated + mStiffnessAlpha * mStiffness;
        // NOTE: Regardless of using Ktime to calculate K, we still need to filter K.
    }

    if (mActuatorCounter == number_of_actuators-1) {
        mActuatorCounter = 0;
    } else{
        mActuatorCounter++;
    }
    KRATOS_CATCH("");
}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::AddTableToSubModelPart(const unsigned int TableId,
                                                                           const Parameters TableParameters,
                                                                           ModelPart* pSubModelPart) {
    KRATOS_TRY
    TableType::Pointer p_table = Kratos::make_shared<TableType>();
    for (IndexType i = 0; i < TableParameters["data"].size(); ++i) {
        p_table->PushBack(TableParameters["data"][i][0].GetDouble(),
                        TableParameters["data"][i][1].GetDouble());
    }
    pSubModelPart->AddTable(TableId, p_table);
    KRATOS_CATCH("");
}

//***************************************************************************************************************

Parameters MultiaxialControlModuleGeneralized2DUtilities::GetDefaultParametersForRadialActuator() {
    KRATOS_TRY
    Parameters default_parameters( R"(
    {
        "Parameters"    : {
            "actuator_name": "Radial",
            "initial_velocity" : 0.0,
            "compression_length" : 2.0,
            "young_modulus" : 7.0e9
        },
        "list_of_dem_boundaries": [],
        "list_of_fem_boundaries": [{
            "model_part_name" : "1",
            "outer_normal": [0.0,0.0,0.0]
        }],
        "target_stress_table": {
            "input_variable": "TIME",
            "output_variable": "TARGET_STRESS",
            "data": [
                [0.0, 0.0],
                [0.7, -1.0e6]
            ]
        }
    }  )" );

    return default_parameters;
    KRATOS_CATCH("");
}

Parameters MultiaxialControlModuleGeneralized2DUtilities::GetDefaultParametersForRadialMultiDofsActuator() {
    KRATOS_TRY
    Parameters default_parameters( R"(
    {
        "Parameters"    : {
            "actuator_name": "Radial",
            "initial_velocity" : 0.0,
            "compression_length" : 2.0,
            "young_modulus" : 7.0e9
        },
        "list_of_dem_boundaries": [],
        "list_of_fem_boundaries": [{
            "model_part_name" : "1",
            "outer_normal": [0.0,0.0,0.0]
        }]
    }  )" );

    return default_parameters;
    KRATOS_CATCH("");
}

Parameters MultiaxialControlModuleGeneralized2DUtilities::GetDefaultParametersForZActuator() {
    KRATOS_TRY
    Parameters default_parameters( R"(
    {
        "Parameters"    : {
            "actuator_name": "Z",
            "initial_velocity" : 0.0,
            "compression_length" : 1.0,
            "young_modulus" : 7.0e9
        },
        "list_of_dem_boundaries": [{
            "model_part_name" : "PartsCont_solid",
            "outer_normal": [0.0,0.0,1.0]
        }],
        "target_stress_table": {
            "input_variable": "TIME",
            "output_variable": "TARGET_STRESS",
            "data": [
                [0.0, 0.0],
                [0.7, -1.0e6]
            ]
        }
    }  )" );

    return default_parameters;
    KRATOS_CATCH("");
}

Parameters MultiaxialControlModuleGeneralized2DUtilities::GetDefaultParametersForXOrYActuator() {
    KRATOS_TRY
    Parameters default_parameters( R"(
    {
        "Parameters"    : {
            "actuator_name": "X",
            "initial_velocity" : 0.0,
            "compression_length" : 0.1524,
            "young_modulus" : 7.0e9
        },
        "list_of_dem_boundaries": [],
        "list_of_fem_boundaries": [{
            "model_part_name" : "left",
            "outer_normal": [-1.0,0.0,0.0]
            },{
            "model_part_name" : "right",
            "outer_normal": [1.0,0.0,0.0]
        }],
        "target_stress_table": {
            "input_variable": "TIME",
            "output_variable": "TARGET_STRESS",
            "data": [
                [0.0, 0.0],
                [5.0e-7, -5.0e4]
            ]
        }
    }  )" );

    return default_parameters;
    KRATOS_CATCH("");
}


void MultiaxialControlModuleGeneralized2DUtilities::CalculateMaximumInputStressVariationRate() {
    KRATOS_TRY
    TableType::Pointer pTargetStressTable;

    mMaximumInputStressVariationRate = -std::numeric_limits<double>::infinity();

    // Loop through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        if (actuator_name == "RadialMultiDofs") { continue; } // RadialMultiDofs target stress is fed externally
        const std::vector<ModelPart*>& FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        const std::vector<ModelPart*>& DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
        unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
        if (actuator_name == "Z") {
            pTargetStressTable = (*(DEMSubModelPartList[0])).pGetTable(target_stress_table_id);
        } else {
            pTargetStressTable = (*(FEMSubModelPartList[0])).pGetTable(target_stress_table_id);
        }
        const auto table_data = pTargetStressTable->Data();

        // Maximum reaction variation rate
        KRATOS_ERROR_IF(table_data.size() < 2) << "Actuator " << actuator_name << " has a target stress table with less than 2 rows." << std::endl;
        double input_variation_rate = mMaximumInputStressVariationRate;
        // Loop through the rows of the target stress table
        for(unsigned int i = 1 ; i < table_data.size() ; i++) {
            if(std::abs(table_data[i].first - table_data[i-1].first) > std::numeric_limits<double>::epsilon()) {
                input_variation_rate = std::abs((table_data[i].second[0] - table_data[i-1].second[0]) / (table_data[i].first - table_data[i-1].first));
            } else {
                KRATOS_ERROR << "Actuator " << actuator_name << " has a target stress table with repeated abscissa input." << std::endl;
            }

            if(input_variation_rate > mMaximumInputStressVariationRate) {
                mMaximumInputStressVariationRate = input_variation_rate;
            }
        }
    }
    KRATOS_CATCH("");
}

void MultiaxialControlModuleGeneralized2DUtilities::SetProvidedInitialVelocityToTheControlledBoundaries() {
    KRATOS_TRY
    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string& actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*>& SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "Radial" || actuator_name == "RadialMultiDofs") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<number_of_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                noalias(r_displacement) = ZeroVector(3);
                noalias(r_delta_displacement) = ZeroVector(3);
                r_velocity[0] = mVelocity[ind] * cos_theta;
                r_velocity[1] = mVelocity[ind] * sin_theta;
                r_velocity[2] = 0.0;
            }
        } else if (actuator_name == "Z") {
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] = 0.0;
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int number_of_nodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<number_of_nodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_displacement) = ZeroVector(3);
                    noalias(r_delta_displacement) = ZeroVector(3);
                    noalias(r_velocity) += mVelocity[ind] * mFEMOuterNormals[actuator_name][i]; //MAC: Additive for corners, IP?
                    r_velocity[2] = 0.0;
                }
            }
        }
    }
    KRATOS_CATCH("");
}
}  // namespace Kratos
