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
#include "custom_utilities/multiaxial_control_module_fem_dem_generalized_2d_utilities.hpp"

namespace Kratos
{

// Before FEM and DEM solution
void MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteInitialize() {
    KRATOS_TRY;

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*> SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        if (actuator_name == "Radial") {
            // In axisymmetric cases we assume there is only 1 actuator in the FEM boundary
            ModelPart& rSubModelPart = *(SubModelPartList[0]);
            // Iterate through nodes of Fem boundary
            const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
            ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
            #pragma omp parallel for
            for(int i = 0; i<NNodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Fix(DISPLACEMENT_X);
                it->Fix(DISPLACEMENT_Y);
                it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
            }
        } else if (actuator_name == "Z") {
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] = 0.0;
        } else if (actuator_name == "X") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->Fix(DISPLACEMENT_X);
                    it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                }
            }
        } else if (actuator_name == "Y") {
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                // Iterate through nodes of Fem boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->Fix(DISPLACEMENT_Y);
                    it->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

//***************************************************************************************************************

// Before FEM and DEM solution
void MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteInitializeSolutionStep() {
    KRATOS_TRY;

    const double current_time = mrFemModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrFemModelPart.GetProcessInfo()[DELTA_TIME];
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
            const std::string actuator_name = mVectorOfActuatorNames[ind];
            std::vector<ModelPart*> FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
            unsigned int target_stress_table_id = mTargetStressTableIds[actuator_name];
            TableType::Pointer pFEMTargetStressTable = (*(FEMSubModelPartList[0])).pGetTable(target_stress_table_id);
            next_target_stress[ind] = pFEMTargetStressTable->GetValue(mCMTime);
        }

        Vector target_stress_perturbation(number_of_actuators);
        noalias(target_stress_perturbation) = this->GetPerturbations(next_target_stress,mCMTime);
        noalias(next_target_stress) += target_stress_perturbation;

        // Calculate velocity
        this->CalculateVelocity(next_target_stress, current_time);
    }

    // Move Actuators

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*> SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
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
                r_displacement[0] += mVelocity[ind] * cos_theta * delta_time;
                r_displacement[1] += mVelocity[ind] * sin_theta * delta_time;
            }
        } else if (actuator_name == "Z") {
            // DEM
            mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE] += mVelocity[ind]*delta_time/1.0;
            // FEM
            const double imposed_z_strain = mrDemModelPart.GetProcessInfo()[IMPOSED_Z_STRAIN_VALUE];
            const ProcessInfo& CurrentProcessInfo = mrFemModelPart.GetProcessInfo();
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < SubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(SubModelPartList[i]);
                int NElems = static_cast<int>(rSubModelPart.Elements().size());
                ModelPart::ElementsContainerType::iterator elem_begin = rSubModelPart.ElementsBegin();
                #pragma omp parallel for
                for(int j = 0; j < NElems; j++)
                {
                    ModelPart::ElementsContainerType::iterator itElem = elem_begin + j;
                    Element::GeometryType& rGeom = itElem->GetGeometry();
                    GeometryData::IntegrationMethod MyIntegrationMethod = itElem->GetIntegrationMethod();
                    const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                    unsigned int NumGPoints = IntegrationPoints.size();
                    std::vector<double> imposed_z_strain_vector(NumGPoints);
                    // Loop through GaussPoints
                    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                    {
                        imposed_z_strain_vector[GPoint] = imposed_z_strain;
                    }
                    itElem->SetValuesOnIntegrationPoints( IMPOSED_Z_STRAIN_VALUE, imposed_z_strain_vector, CurrentProcessInfo );
                }
            }
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
                    noalias(r_displacement) += mVelocity[ind] * mFEMOuterNormals[actuator_name][i] * delta_time;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

//***************************************************************************************************************

// After FEM and DEM solution
void MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteFinalizeSolutionStep() {
    const double current_time = mrFemModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrFemModelPart.GetProcessInfo()[DELTA_TIME];
    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();

    // Update ReactionStresses
    Vector reaction_stress_estimated(number_of_actuators);
    noalias(reaction_stress_estimated) = this->MeasureReactionStress(REACTION);
    noalias(mReactionStress) = (1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress;
    // noalias(mReactionStress) = 1.0/(1.0 - std::pow(mReactionAlpha,mStep)) *
    //                         ((1.0 - mReactionAlpha) * reaction_stress_estimated + mReactionAlpha * mReactionStress);

    noalias(mDisplacement) += mVelocity * mCMDeltaTime;

    // Update Stiffness matrix
    if (current_time > (mKTime - 0.5 * delta_time)) {

        // Update K if DeltaDisplacement is invertible
        this->CalculateStiffness();
    }

    // Print results

    // Iterate through all on plane actuators to set variables to 0.0
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*> SubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
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
                    array_1d<double,3>& r_target_stress = it->GetValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->GetValue(REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->GetValue(LOADING_VELOCITY);
                    noalias(r_target_stress) = ZeroVector(3);
                    noalias(r_reaction_stress) = ZeroVector(3);
                    noalias(r_loading_velocity) = ZeroVector(3);
                }
            }
        }
    }

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*> FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        std::vector<ModelPart*> DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
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
                it->GetValue(TARGET_STRESS_X) = current_target_stress * cos_theta;
                it->GetValue(TARGET_STRESS_Y) = current_target_stress * sin_theta;
                it->GetValue(REACTION_STRESS_X) = mReactionStress[ind] * cos_theta;
                it->GetValue(REACTION_STRESS_Y) = mReactionStress[ind] * sin_theta;
                it->GetValue(LOADING_VELOCITY_X) = mVelocity[ind] * cos_theta;
                it->GetValue(LOADING_VELOCITY_Y) = mVelocity[ind] * sin_theta;
            }
        } else if (actuator_name == "Z") {
            // Note: we only print on FEM nodes
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                // Iterate through nodes of FEM boundary
                const int NNodes = static_cast<int>(rSubModelPart.Nodes().size());
                ModelPart::NodesContainerType::iterator it_begin = rSubModelPart.NodesBegin();
                TableType::Pointer TargetStressTable = rSubModelPart.pGetTable(target_stress_table_id);
                double current_target_stress = TargetStressTable->GetValue(current_time);
                #pragma omp parallel for
                for(int j = 0; j<NNodes; j++) {
                    ModelPart::NodesContainerType::iterator it = it_begin + j;
                    it->GetValue(TARGET_STRESS_Z) = current_target_stress;
                    it->GetValue(REACTION_STRESS_Z) = mReactionStress[ind];
                    it->GetValue(LOADING_VELOCITY_Z) = mVelocity[ind];
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
                    array_1d<double,3>& r_target_stress = it->GetValue(TARGET_STRESS);
                    array_1d<double,3>& r_reaction_stress = it->GetValue(REACTION_STRESS);
                    array_1d<double,3>& r_loading_velocity = it->GetValue(LOADING_VELOCITY);
                    noalias(r_target_stress) += current_target_stress * mFEMOuterNormals[actuator_name][i];
                    noalias(r_reaction_stress) += mReactionStress[ind] * mFEMOuterNormals[actuator_name][i];
                    noalias(r_loading_velocity) += mVelocity[ind] * mFEMOuterNormals[actuator_name][i];
                }
            }
        }
    }
}

//***************************************************************************************************************

Vector MultiaxialControlModuleFEMDEMGeneralized2DUtilities::MeasureReactionStress(const Variable<array_1d<double,3>>& rVariable) {

    const unsigned int number_of_actuators = mListsOfFEMSubModelPartsForEachActuator.size();
    Vector reaction_stress(number_of_actuators);
    noalias(reaction_stress) = ZeroVector(number_of_actuators);

    // Iterate through all actuators
    for(unsigned int ind = 0; ind < mVectorOfActuatorNames.size(); ind++) {
        const std::string actuator_name = mVectorOfActuatorNames[ind];
        std::vector<ModelPart*> FEMSubModelPartList = mListsOfFEMSubModelPartsForEachActuator[actuator_name];
        std::vector<ModelPart*> DEMSubModelPartList = mListsOfDEMSubModelPartsForEachActuator[actuator_name];
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
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(rVariable);
                    // Unit normal vector pointing outwards
                    array_1d<double,3> radial_normal;
                    radial_normal[0] = it->X();
                    radial_normal[1] = it->Y();
                    radial_normal[2] = 0.0;
                    double inv_norm = 1.0/norm_2(radial_normal);
                    radial_normal[0] *= inv_norm;
                    radial_normal[1] *= inv_norm;
                    face_reaction += inner_prod(r_force,radial_normal);
                }
            }
            if (std::abs(face_area) > 1.0e-12) {
                reaction_stress[ind] = face_reaction/face_area;
            } else {
                reaction_stress[ind] = 0.0;
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
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                int NElems = static_cast<int>(rSubModelPart.Elements().size());
                ModelPart::ElementsContainerType::iterator elem_begin = rSubModelPart.ElementsBegin();
                #pragma omp parallel for reduction(+:face_area)
                for(int j = 0; j < NElems; j++)
                {
                    ModelPart::ElementsContainerType::iterator itElem = elem_begin + j;
                    face_area += itElem->GetGeometry().Area();
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
            const ProcessInfo& CurrentProcessInfo = mrFemModelPart.GetProcessInfo();
            // Iterate through all FEMBoundaries
            for (unsigned int i = 0; i < FEMSubModelPartList.size(); i++) {
                ModelPart& rSubModelPart = *(FEMSubModelPartList[i]);
                int NElems = static_cast<int>(rSubModelPart.Elements().size());
                ModelPart::ElementsContainerType::iterator elem_begin = rSubModelPart.ElementsBegin();
                #pragma omp parallel for reduction(+:face_reaction)
                for(int j = 0; j < NElems; j++)
                {
                    ModelPart::ElementsContainerType::iterator itElem = elem_begin + j;
                    Element::GeometryType& rGeom = itElem->GetGeometry();
                    GeometryData::IntegrationMethod MyIntegrationMethod = itElem->GetIntegrationMethod();
                    const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                    unsigned int NumGPoints = IntegrationPoints.size();
                    std::vector<Vector> stress_vector(NumGPoints);
                    itElem->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, CurrentProcessInfo );
                    const double area_over_gp = rGeom.Area()/NumGPoints;
                    // Loop through GaussPoints
                    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                    {
                        face_reaction += stress_vector[GPoint][2] * area_over_gp;
                    }
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
                    array_1d<double,3>& r_force = it->FastGetSolutionStepValue(rVariable);
                    face_reaction += inner_prod(r_force,mFEMOuterNormals[actuator_name][i]);
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
}

}  // namespace Kratos
