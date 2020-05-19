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

// Before DEM solution
void MultiaxialControlModuleGeneralized2DUtilities::ExecuteInitialize() {
    KRATOS_TRY;

    // Iterate through all on plane actuators to set velocities to 0.0
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> SubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        if (actuator_name != "Radial" && actuator_name != "Z") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_velocity) = ZeroVector(3);
                }
            }
        }
    }

    // Iterate through all actuators
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> SubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                noalias(r_displacement) = ZeroVector(3);
                noalias(r_delta_displacement) = ZeroVector(3);
                r_velocity[0] = mVelocity[map_index] * cos_theta;
                r_velocity[1] = mVelocity[map_index] * sin_theta;
                r_velocity[2] = 0.0;
            }
        } else if (actuator_name == "Z") {
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] = 0.0;
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_displacement) = ZeroVector(3);
                    noalias(r_delta_displacement) = ZeroVector(3);
                    noalias(r_velocity) += mVelocity[map_index] * mFEMOuterNormals[actuator_name][i];
                    r_velocity[2] = 0.0;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

//***************************************************************************************************************

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

        const unsigned int number_of_actuators = mFEMBoundariesSubModelParts.size();

        Vector next_target_stress(number_of_actuators);
        noalias(next_target_stress) = ZeroVector(number_of_actuators);

        // Iterate through all actuators
        for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
            const std::string actuator_name = mOrderedMapKeys[map_index];
            std::vector<ModelPart*> FEMSubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
            std::vector<ModelPart*> DEMSubModelPartList = mDEMBoundariesSubModelParts[actuator_name];
            unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
            if (actuator_name == "Z") {
                TableType::Pointer pDEMTargetStressTable = (*(DEMSubModelPartList[0])).pGetTable(target_stress_table_id);
                next_target_stress[map_index] = pDEMTargetStressTable->GetValue(mCMTime);
            } else {
                TableType::Pointer pFEMTargetStressTable = (*(FEMSubModelPartList[0])).pGetTable(target_stress_table_id);
                next_target_stress[map_index] = pFEMTargetStressTable->GetValue(mCMTime);
            }
        }

        Vector target_stress_perturbation(number_of_actuators);
        noalias(target_stress_perturbation) = GetPerturbations(next_target_stress,mCMTime);
        noalias(next_target_stress) += target_stress_perturbation;

        // Calculate velocity
        CalculateVelocity(next_target_stress, current_time);
    }

    // Move Actuators

    // Iterate through all actuators
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> SubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                r_velocity[0] = mVelocity[map_index] * cos_theta;
                r_velocity[1] = mVelocity[map_index] * sin_theta;
                r_velocity[2] = 0.0;
                noalias(r_delta_displacement) = r_velocity * delta_time;
                noalias(r_displacement) += r_delta_displacement;
                noalias(it->Coordinates()) = it->GetInitialPosition().Coordinates() + r_displacement;
            }
        } else if (actuator_name == "Z") {
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] += mVelocity[map_index]*delta_time/1.0;
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& r_delta_displacement = it->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                    array_1d<double,3>& r_velocity = it->FastGetSolutionStepValue(VELOCITY);
                    noalias(r_velocity) = mVelocity[map_index] * mFEMOuterNormals[actuator_name][i];
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

//***************************************************************************************************************

// After DEM solution
void MultiaxialControlModuleGeneralized2DUtilities::ExecuteFinalizeSolutionStep() {
    const double current_time = mrDemModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrDemModelPart.GetProcessInfo()[DELTA_TIME];
    const unsigned int number_of_actuators = mFEMBoundariesSubModelParts.size();

    // Update ReactionStresses
    Vector reaction_stress_estimated(number_of_actuators);
    noalias(reaction_stress_estimated) = MeasureReactionStress();
    noalias(mReactionStress) = (1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress;
    // noalias(mReactionStress) = 1.0/(1.0 - std::pow(mReactionAlpha,mStep)) * 
    //                         ((1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress);

    // Update Stiffness matrix
    if (current_time > (mCMTime - 0.5 * delta_time)) {

        // Update K if DeltaDisplacement is invertible
        CalculateStiffness();
    }

    // Print results

    // Iterate through all on plane actuators to set variables to 0.0
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> SubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        if (actuator_name != "Radial" && actuator_name != "Z") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_target_stress = it->FastGetSolutionStepValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->FastGetSolutionStepValue(REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->FastGetSolutionStepValue(LOADING_VELOCITY);
                    noalias(r_target_stress) = ZeroVector(3);
                    noalias(r_reaction_stress) = ZeroVector(3);
                    noalias(r_loading_velocity) = ZeroVector(3);
                }
            }
        }
    }

    // Iterate through all actuators
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> FEMSubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        std::vector<ModelPart*> DEMSubModelPartList = mDEMBoundariesSubModelParts[actuator_name];
        unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(FEMSubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
            double current_target_stress = TargetStressTable->GetValue(current_time);
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                const double external_radius = std::sqrt(it->X()*it->X() + it->Y()*it->Y());
                const double cos_theta = it->X()/external_radius;
                const double sin_theta = it->Y()/external_radius;
                it->FastGetSolutionStepValue(TARGET_STRESS_X) = current_target_stress * cos_theta;
                it->FastGetSolutionStepValue(TARGET_STRESS_Y) = current_target_stress * sin_theta;
                it->FastGetSolutionStepValue(REACTION_STRESS_X) = mReactionStress[map_index] * cos_theta;
                it->FastGetSolutionStepValue(REACTION_STRESS_Y) = mReactionStress[map_index] * sin_theta;
                it->FastGetSolutionStepValue(LOADING_VELOCITY_X) = mVelocity[map_index] * cos_theta;
                it->FastGetSolutionStepValue(LOADING_VELOCITY_Y) = mVelocity[map_index] * sin_theta;
            }
        } else if (actuator_name == "Z") {
            // Iterate through all DEMBoundaries
            for (unsigned int i = 0; i < DEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(DEMSubModelPartList[i]);
                // Iterate through nodes of DEM boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
                double current_target_stress = TargetStressTable->GetValue(current_time);
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->FastGetSolutionStepValue(TARGET_STRESS_Z) = current_target_stress;
                    it->FastGetSolutionStepValue(REACTION_STRESS_Z) = mReactionStress[map_index];
                    it->FastGetSolutionStepValue(LOADING_VELOCITY_Z) = mVelocity[map_index];
                }
                mrDemModelPart.GetProcessInfo()[TARGET_STRESS_Z] = std::abs(current_target_stress);
            }
        } else {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
                double current_target_stress = TargetStressTable->GetValue(current_time);
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_target_stress = it->FastGetSolutionStepValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->FastGetSolutionStepValue(REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->FastGetSolutionStepValue(LOADING_VELOCITY);
                    noalias(r_target_stress) += current_target_stress * mFEMOuterNormals[actuator_name][i];
                    noalias(r_reaction_stress) += mReactionStress[map_index] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_loading_velocity) += mVelocity[map_index] * mFEMOuterNormals[actuator_name][i];
                }
            }
        }
    }
}

//***************************************************************************************************************

Vector MultiaxialControlModuleGeneralized2DUtilities::MeasureReactionStress() {

    const unsigned int number_of_actuators = mFEMBoundariesSubModelParts.size();
    Vector reaction_stress(number_of_actuators);
    noalias(reaction_stress) = ZeroVector(number_of_actuators);

    // Iterate through all actuators
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        std::vector<ModelPart*> FEMSubModelPartList = mFEMBoundariesSubModelParts[actuator_name];
        std::vector<ModelPart*> DEMSubModelPartList = mDEMBoundariesSubModelParts[actuator_name];
        double face_area = 0.0;
        double face_reaction = 0.0;
        if (actuator_name == "Radial") {
            // Calculate face_area
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through conditions of FEM boundary
                const int NCons = static_cast<int>(rSubModelPart.Conditions().size());
                ModelPart::ConditionsContainerType::iterator con_begin = rSubModelPart.ConditionsBegin();
                #pragma omp parallel for reduction(+:face_area)
                for(int j = 0; j < NCons; j++) {
                    ModelPart::ConditionsContainerType::iterator itCond = con_begin + j;
                    face_area += itCond->GetGeometry().Area();
                }
            }
            // Calculate face_reaction
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for reduction(+:face_reaction)
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(CONTACT_FORCES);
                    // Unit normal vector pointing outwards
                    array_1d<double,3> radial_normal;
                    radial_normal[0] = it->X();
                    radial_normal[1] = it->Y();
                    radial_normal[2] = 0.0;
                    double inv_norm = 1.0/norm_2(radial_normal);
                    radial_normal[0] *= inv_norm;
                    radial_normal[1] *= inv_norm;
                    face_reaction -= inner_prod(r_force,radial_normal);
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[map_index] = face_reaction/face_area;
            } else {
                reaction_stress[map_index] = 0.0;
            }
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
                reaction_stress[map_index] = face_reaction/face_area;
            } else {
                reaction_stress[map_index] = 0.0;
            }
        } else {
            // Calculate face_area
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through conditions of FEM boundary
                const int NCons = static_cast<int>(rSubModelPart.Conditions().size());
                ModelPart::ConditionsContainerType::iterator con_begin = rSubModelPart.ConditionsBegin();
                #pragma omp parallel for reduction(+:face_area)
                for(int j = 0; j < NCons; j++) {
                    ModelPart::ConditionsContainerType::iterator itCond = con_begin + j;
                    face_area += itCond->GetGeometry().Area();
                }
            }
            // Calculate face_reaction
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for reduction(+:face_reaction)
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(CONTACT_FORCES);
                    face_reaction -= inner_prod(r_force,mFEMOuterNormals[actuator_name][i]);
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[map_index] = face_reaction/face_area;
            } else {
                reaction_stress[map_index] = 0.0;
            }
        }
    }

    return reaction_stress;
}

//***************************************************************************************************************

Vector MultiaxialControlModuleGeneralized2DUtilities::GetPerturbations(const Vector& rTargetStress, const double& rTime) {

    const unsigned int number_of_actuators = rTargetStress.size();
    Vector stress_perturbation(number_of_actuators);
    noalias(stress_perturbation) = ZeroVector(number_of_actuators);

    // Iterate through all actuators
    for(unsigned int map_index = 0; map_index < mOrderedMapKeys.size(); map_index++) {
        const std::string actuator_name = mOrderedMapKeys[map_index];
        if (actuator_name == "Z") {
            stress_perturbation[map_index] = 0.0;
        } else {
            double amplitude = rTargetStress[map_index] * mPerturbationTolerance;
            double omega = 2.0 * Globals::Pi / (mPerturbationPeriod * mCMDeltaTime);
            double phi = map_index * 2.0 * Globals::Pi / number_of_actuators;
            stress_perturbation[map_index] = amplitude * std::sin(omega * rTime + phi);
        }
    }
    return stress_perturbation;
}

//***************************************************************************************************************

double MultiaxialControlModuleGeneralized2DUtilities::GetConditionNumber(const Matrix& rInputMatrix, const Matrix& rInvertedMatrix) {

    // Find the condition number to define is inverse is OK
    const double input_matrix_norm = norm_frobenius(rInputMatrix);
    const double inverted_matrix_norm = norm_frobenius(rInvertedMatrix);

    const double cond_number = input_matrix_norm * inverted_matrix_norm ;

    return cond_number;
}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::CalculateVelocity(const Vector& r_next_target_stress, const double& r_current_time) {

    const unsigned int number_of_actuators = mFEMBoundariesSubModelParts.size();

    Vector delta_target_stress(number_of_actuators);
    noalias(delta_target_stress) = r_next_target_stress-mReactionStress;
    Vector velocity_estimated(number_of_actuators);

    if (mMultiAxial == true) {
        Matrix k_inverse(number_of_actuators,number_of_actuators);
        double k_det = 0.0;
        MathUtils<double>::InvertMatrix(mStiffness, k_inverse, k_det, -1.0);
        // TODO: check tolerance in CheckConditionNumber
        const bool is_k_invertible = MathUtils<double>::CheckConditionNumber(mStiffness, 
                                                        k_inverse, std::numeric_limits<double>::epsilon(),
                                                        false);
        const double k_condition_number = GetConditionNumber(mStiffness,k_inverse);

        KRATOS_WATCH("Begin Updating velocity.........")
        KRATOS_WATCH(r_current_time)
        KRATOS_WATCH(mCMTime)
        KRATOS_WATCH(mStiffness)
        KRATOS_WATCH(k_condition_number)
        KRATOS_WATCH(mVelocity)

        Vector velocity_perturbation(number_of_actuators);
        noalias(velocity_perturbation) = GetPerturbations(mLimitVelocities,r_current_time);
        if (is_k_invertible == false || std::isnan(k_condition_number)) {
            noalias(velocity_estimated) = mVelocity + velocity_perturbation;
            std::cout << "Stiffness matrix is not invertible. Keeping loading velocity constant" << std::endl;            
        } else {
            noalias(velocity_estimated) = prod(k_inverse,delta_target_stress)/mCMDeltaTime;
        }
        noalias(mVelocity) = (1.0 - mVelocityAlpha) * velocity_estimated + mVelocityAlpha * mVelocity;

        KRATOS_WATCH("Updating velocity.........")
        KRATOS_WATCH(mVelocity)

        for (unsigned int i = 0; i < mVelocity.size(); i++) {
            if (std::abs(mVelocity[i]) > std::abs(mLimitVelocities[i])) {
                if (mVelocity[i] > 0.0) {
                    mVelocity[i] = std::abs(mLimitVelocities[i]) + std::abs(velocity_perturbation[i]);
                } else {
                    mVelocity[i] = - std::abs(mLimitVelocities[i]) - std::abs(velocity_perturbation[i]);
                }
            }
        }

        KRATOS_WATCH("End Updating velocity.........")
        KRATOS_WATCH(mVelocity)
    } else {
        for (unsigned int i = 0; i < mVelocity.size(); i++) {
            velocity_estimated[i] = delta_target_stress[i]/(mStiffness(i,i)*mCMDeltaTime);
        }
        noalias(mVelocity) = (1.0 - mVelocityAlpha) * velocity_estimated + mVelocityAlpha * mVelocity;

        for (unsigned int i = 0; i < mVelocity.size(); i++) {
            if (std::abs(mVelocity[i]) > std::abs(mLimitVelocities[i])) {
                if (mVelocity[i] > 0.0) {
                    mVelocity[i] = std::abs(mLimitVelocities[i]);
                } else {
                    mVelocity[i] = - std::abs(mLimitVelocities[i]);
                }
            }
        }
    }


}

//***************************************************************************************************************

void MultiaxialControlModuleGeneralized2DUtilities::CalculateStiffness() {

    const unsigned int number_of_actuators = mFEMBoundariesSubModelParts.size();

    Vector delta_reaction_stress(number_of_actuators);
    noalias(delta_reaction_stress) = mReactionStress - mReactionStressOld;
    noalias(mReactionStressOld) = mReactionStress;

    for (unsigned int i = 0; i < number_of_actuators; i++) {
        mDeltaDisplacement(i,mActuatorCounter) = mVelocity[i]*mCMDeltaTime;
        mDeltaReactionStress(i,mActuatorCounter) = delta_reaction_stress[i];
    }

    Matrix k_estimated(number_of_actuators,number_of_actuators);

    if (mMultiAxial == true) {
        if (mCMStep > number_of_actuators-1) {
            Matrix delta_displacement_inverse(number_of_actuators,number_of_actuators);
            double delta_displacement_det = 0.0;
            MathUtils<double>::InvertMatrix(mDeltaDisplacement, delta_displacement_inverse, delta_displacement_det,-1.0);
            // TODO: check tolerance in CheckConditionNumber
            const bool is_delta_displacement_invertible = MathUtils<double>::CheckConditionNumber(mDeltaDisplacement, 
                                                            delta_displacement_inverse, 1.0e-10, 
                                                            false);
            const double delta_displacement_condition_number = GetConditionNumber(mDeltaDisplacement,delta_displacement_inverse);

            KRATOS_WATCH("Begin Updating K.........")
            KRATOS_WATCH(mrDemModelPart.GetProcessInfo()[TIME])
            KRATOS_WATCH(mCMTime)
            KRATOS_WATCH(mDeltaDisplacement)
            KRATOS_WATCH(delta_displacement_condition_number)
            KRATOS_WATCH(mStiffness)


            if (is_delta_displacement_invertible == false || std::isnan(delta_displacement_condition_number) ) {
                noalias(k_estimated) = mStiffness;
                std::cout << "Delta displacement matrix is not invertible. Keeping stiffness matrix constant" << std::endl;
            } else {
                noalias(k_estimated) = prod(mDeltaReactionStress,delta_displacement_inverse);
            }
            noalias(mStiffness) = (1.0 - mStiffnessAlpha) * k_estimated + mStiffnessAlpha * mStiffness;

            KRATOS_WATCH("End Updating K........")
            KRATOS_WATCH(mStiffness)    
        }
    } else {
        for (unsigned int i = 0; i < number_of_actuators; i++) {
            k_estimated(i,i) = mDeltaReactionStress(i,mActuatorCounter)/mDeltaDisplacement(i,mActuatorCounter);
        }
        noalias(mStiffness) = (1.0 - mStiffnessAlpha) * k_estimated + mStiffnessAlpha * mStiffness;
    }

    if (mActuatorCounter == number_of_actuators-1) {
        mActuatorCounter = 0;
    } else{
        mActuatorCounter++;
    }
}


}  // namespace Kratos
