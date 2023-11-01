/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu, chengshun.shang1996@gmail.com
// Date: June 2023
/////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <algorithm>

// Project includes
#include "DEM_smooth_joint_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos{

DEMContinuumConstitutiveLaw::Pointer DEM_smooth_joint::Clone() const{
    DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_smooth_joint(*this));
    return p_clone;
}

void DEM_smooth_joint::Initialize(SphericContinuumParticle* element1,
                                                 SphericContinuumParticle* element2,
                                                 Properties::Pointer pProps) {
        mpProperties = pProps;
        double GlobalJointNormal[3] = {0.0};
        GlobalJointNormal[0] = (*mpProperties)[JOINT_NORMAL_DIRECTION_X];
        GlobalJointNormal[1] = (*mpProperties)[JOINT_NORMAL_DIRECTION_Y];
        GlobalJointNormal[2] = (*mpProperties)[JOINT_NORMAL_DIRECTION_Z];
        array_1d<double, 3> GlobalOtherToMeVector;
        array_1d<double, 3> LocalOtherToMeVector;
        noalias(GlobalOtherToMeVector) = element1->GetGeometry()[0].Coordinates() - element2->GetGeometry()[0].Coordinates();
        double LocalCoordSystem[3][3];
        double Distance = DEM_MODULUS_3(GlobalOtherToMeVector);
        GeometryFunctions::ComputeContactLocalCoordSystem(GlobalOtherToMeVector, Distance, LocalCoordSystem);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalJointNormal, mLocalJointNormal);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalOtherToMeVector, LocalOtherToMeVector);
        double temp_OtherToMeVector[3] = {0.0};
        temp_OtherToMeVector[0] = LocalOtherToMeVector[0];
        temp_OtherToMeVector[1] = LocalOtherToMeVector[1];
        temp_OtherToMeVector[2] = LocalOtherToMeVector[2];
        double temp_InitialDistanceJoint = GeometryFunctions::DotProduct(temp_OtherToMeVector, mLocalJointNormal);
        mInitialDistanceJoint = std::abs(temp_InitialDistanceJoint);
    }

std::string DEM_smooth_joint::GetTypeOfLaw() {
        std::string type_of_law = "smooth_joint_CL";
        return type_of_law;
    }

//*************************************
// Parameters preparation
//*************************************

void DEM_smooth_joint::TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp){
    BaseClassType::TransferParametersToProperties(parameters, pProp);
}

void DEM_smooth_joint::Check(Properties::Pointer pProp) const {
    
    //two parts: discontinuum part and continuum part

    //*********discontinuum part **************

    if(!pProp->Has(STATIC_FRICTION)){
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable STATIC_FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(STATIC_FRICTION) = 0.0;
    }

    if(!pProp->Has(DYNAMIC_FRICTION)){
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable DYNAMIC_FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
    }

    //********** continuum part ***************
    //this version of smooth joint model is based on [DEM_parallel_bond_CL]

    if(!pProp->Has(JOINT_NORMAL_STIFFNESS)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_NORMAL_STIFFNESS should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_NORMAL_STIFFNESS) = 1e9;
    }

    if(!pProp->Has(JOINT_TANGENTIAL_STIFFNESS)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_TANGENTIAL_STIFFNESS should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_TANGENTIAL_STIFFNESS) = 1e9;
    }
    
    if(!pProp->Has(BOND_SIGMA_MAX)) {

        KRATOS_ERROR << "Variable BOND_SIGMA_MAX should be present in the properties when using DEM_smooth_joint_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_SIGMA_MAX_DEVIATION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_SIGMA_MAX_DEVIATION should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_SIGMA_MAX_DEVIATION) = 0.0;
    }

    if(!pProp->Has(BOND_TAU_ZERO)) {

        KRATOS_ERROR << "Variable BOND_TAU_ZERO should be present in the properties when using DEM_smooth_joint_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_TAU_ZERO_DEVIATION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_TAU_ZERO_DEVIATION should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_TAU_ZERO_DEVIATION) = 0.0;
    }

    if(!pProp->Has(BOND_INTERNAL_FRICC)) {

        KRATOS_ERROR << "Variable BOND_INTERNAL_FRICC should be present in the properties when using DEM_smooth_joint_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_RADIUS_FACTOR)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_RADIUS_FACTOR should be present in the properties when using DEM_smooth_joint_CL. 1.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_RADIUS_FACTOR) = 1.0;
    }

    if(!pProp->Has(JOINT_NORMAL_DIRECTION_X)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_NORMAL_DIRECTION_X should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_NORMAL_DIRECTION_X) = 0.0;
    }

    if(!pProp->Has(JOINT_NORMAL_DIRECTION_Y)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_NORMAL_DIRECTION_Y should be present in the properties when using DEM_smooth_joint_CL. 1.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_NORMAL_DIRECTION_Y) = 1.0;
    }

    if(!pProp->Has(JOINT_NORMAL_DIRECTION_Z)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_NORMAL_DIRECTION_Z should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_NORMAL_DIRECTION_Z) = 0.0;
    }

    if(!pProp->Has(JOINT_FRICTION_COEFF)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable JOINT_FRICTION_COEFF should be present in the properties when using DEM_smooth_joint_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(JOINT_FRICTION_COEFF) = 0.0;
    }

    if (!pProp->Has(IS_UNBREAKABLE)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable IS_UNBREAKABLE was not present in the properties when using DEM_smooth_joint_CL. False value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(IS_UNBREAKABLE) = false;
    }
} // CHECK()

// TODO: Calculate Contact Area = Calculate Bond area? --> Yes
void DEM_smooth_joint::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

    KRATOS_TRY

    const double bond_radius_factor = (*mpProperties)[BOND_RADIUS_FACTOR];
    double bond_radius = std::min(radius, other_radius) * bond_radius_factor;
    calculation_area = Globals::Pi * bond_radius * bond_radius;

    KRATOS_CATCH("")
}

double DEM_smooth_joint::CalculateContactArea(double radius, double other_radius, Vector& v) {
    
    KRATOS_TRY

    double a = 0.0;
    CalculateContactArea(radius, other_radius, a);
    unsigned int old_size = v.size();
    Vector backup = v;
    v.resize(old_size + 1, false);
    v[old_size] = a;
    for (unsigned int i=0; i<old_size; i++) {
        v[i] = backup[i];
    }
    return a;

    KRATOS_CATCH("")
}

void DEM_smooth_joint::GetContactArea(const double radius, const double other_radius, const Vector& vector_of_initial_areas, const int neighbour_position, double& calculation_area) {
    
    KRATOS_TRY

    CalculateContactArea(radius, other_radius, calculation_area);

    KRATOS_CATCH("")
}


// TODO: In this function, it is better to replace 'kn_el' with 'kn_bond' and 'kt_el' with 'kt_bond'.
void DEM_smooth_joint::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

    KRATOS_TRY

    //for bonded part
    //const double bond_equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];
    //kn_el = bond_equiv_young * calculation_area / mInitialDistanceJoint;
    //kt_el = kn_el / (*mpProperties)[BOND_KNKS_RATIO];
    kn_el = (*mpProperties)[JOINT_NORMAL_STIFFNESS];
    kt_el = (*mpProperties)[JOINT_TANGENTIAL_STIFFNESS];

    double GlobalJointNormal[3] = {0.0};
    GlobalJointNormal[0] = (*mpProperties)[JOINT_NORMAL_DIRECTION_X]; // X = -1 * sin(alpha) 
    GlobalJointNormal[1] = (*mpProperties)[JOINT_NORMAL_DIRECTION_Y]; // Y = cos(alpha)
    GlobalJointNormal[2] = (*mpProperties)[JOINT_NORMAL_DIRECTION_Z]; // Z = 0
    array_1d<double, 3> OtherToMeVector;
    noalias(OtherToMeVector) = element1->GetGeometry()[0].Coordinates() - element2->GetGeometry()[0].Coordinates();
    double LocalCoordSystem[3][3];
    double Distance = DEM_MODULUS_3(OtherToMeVector);
    GeometryFunctions::ComputeContactLocalCoordSystem(OtherToMeVector, Distance, LocalCoordSystem);
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalJointNormal, mLocalJointNormal);

    KRATOS_CATCH("")

}

double DEM_smooth_joint::LocalMaxSearchDistance(const int i,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2) {
    KRATOS_TRY

    double tension_limit;

    // calculation of equivalent Young modulus
    //const double& equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];

    const double my_radius = element1->GetRadius();
    const double other_radius = element2->GetRadius();
    double calculation_area = 0.0;

    Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);
    GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);

    double radius_sum = my_radius + other_radius;
    //double initial_delta = element1->GetInitialDelta(i);
    //double initial_dist = radius_sum - initial_delta;

    // calculation of elastic constants
    //double kn_el = equiv_young * calculation_area / mInitialDistanceJoint;
    double kn_el = (*mpProperties)[JOINT_NORMAL_STIFFNESS]; 

    //tension_limit = GetContactSigmaMax();
    tension_limit = (*mpProperties)[BOND_SIGMA_MAX]; //TODO: add BOND_SIGMA_MAX_DEVIATION

    const double max_normal_bond_force = tension_limit * calculation_area;
    double u_max = max_normal_bond_force / kn_el;

    //TODO: need to choose whether the [if] below is necessary
    if (u_max > 2.0 * radius_sum) {
        u_max = 2.0 * radius_sum;
    } // avoid error in special cases with too high tensile
    

    return u_max;

    KRATOS_CATCH("")
}

//*************************************
// Force calculation
//*************************************

void DEM_smooth_joint::CalculateForces(const ProcessInfo& r_process_info,
                            double OldLocalElasticContactForce[3],
                            double LocalElasticContactForce[3],
                            double LocalElasticExtraContactForce[3],
                            double LocalCoordSystem[3][3],
                            double LocalDeltDisp[3],
                            const double kn_el,
                            const double kt_el,
                            double& contact_sigma,
                            double& contact_tau,
                            double& failure_criterion_state,
                            double equiv_young,
                            double equiv_shear,
                            double indentation,
                            double calculation_area,
                            double& acumulated_damage,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2,
                            int i_neighbour_count,
                            int time_steps,
                            bool& sliding,
                            double &equiv_visco_damp_coeff_normal,
                            double &equiv_visco_damp_coeff_tangential,
                            double LocalRelVel[3],
                            double ViscoDampingLocalContactForce[3]) {
                           
    KRATOS_TRY

    CalculateNormalForces(LocalElasticContactForce,
                        kn_el,
                        equiv_young,
                        indentation,
                        calculation_area,
                        acumulated_damage,
                        element1,
                        element2,
                        i_neighbour_count,
                        time_steps,
                        r_process_info,
                        contact_sigma);

    // Tangential forces are calculated after the viscodamping because the frictional limit bounds the sum of elastic plus viscous forces
    CalculateTangentialForces(OldLocalElasticContactForce,
                            LocalElasticContactForce,
                            LocalElasticExtraContactForce,
                            ViscoDampingLocalContactForce,
                            LocalCoordSystem,
                            LocalDeltDisp,
                            LocalRelVel,
                            kt_el,
                            equiv_shear,
                            contact_tau,
                            indentation,
                            calculation_area,
                            failure_criterion_state,
                            element1,
                            element2,
                            i_neighbour_count,
                            sliding,
                            r_process_info,
                            time_steps);

    KRATOS_CATCH("") 
}

void DEM_smooth_joint::CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps,
                const ProcessInfo& r_process_info,
                double& contact_sigma) {

    KRATOS_TRY
    
    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    array_1d<double, 3> GlobalOtherToMeVector;
    array_1d<double, 3> LocalOtherToMeVector;
    noalias(GlobalOtherToMeVector) = element1->GetGeometry()[0].Coordinates() - element2->GetGeometry()[0].Coordinates();
    double LocalCoordSystem[3][3];
    double Distance = DEM_MODULUS_3(GlobalOtherToMeVector);
    GeometryFunctions::ComputeContactLocalCoordSystem(GlobalOtherToMeVector, Distance, LocalCoordSystem);
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalOtherToMeVector, LocalOtherToMeVector);
    double temp_CurrentOtherToMeVector[3] = {0.0};
    temp_CurrentOtherToMeVector[0] = LocalOtherToMeVector[0];
    temp_CurrentOtherToMeVector[1] = LocalOtherToMeVector[1];
    temp_CurrentOtherToMeVector[2] = LocalOtherToMeVector[2];
    double temp_distance = GeometryFunctions::DotProduct(temp_CurrentOtherToMeVector, mLocalJointNormal);
    const double joint_indentation = mInitialDistanceJoint - std::abs(temp_distance);                                                                                               

    double JointLocalElasticContactForce2 = 0.0;

    if (!failure_type){ //if the bond is not broken
        JointLocalElasticContactForce2 = kn_el * joint_indentation;
    } else { //else the bond is broken
        JointLocalElasticContactForce2 = 0.0;
    }
        
    if(calculation_area){
        contact_sigma = JointLocalElasticContactForce2 / calculation_area;
    }

    LocalElasticContactForce[2] = JointLocalElasticContactForce2;

    KRATOS_CATCH("")  

} // CalculateNormalForces

void DEM_smooth_joint::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                        double LocalElasticContactForce[3],
                        double LocalElasticExtraContactForce[3],
                        double ViscoDampingLocalContactForce[3],
                        double LocalCoordSystem[3][3],
                        double LocalDeltDisp[3],
                        double LocalRelVel[3],
                        const double kt_el,
                        const double equiv_shear,
                        double& contact_tau,
                        double indentation,
                        double calculation_area,
                        double& failure_criterion_state,
                        SphericContinuumParticle* element1,
                        SphericContinuumParticle* element2,
                        int i_neighbour_count,
                        bool& sliding,
                        const ProcessInfo& r_process_info,
                        int time_steps) {

    KRATOS_TRY

    // TODO: add [BOND_TAU_ZERO_DEVIATION]
    //const double& bond_tau_zero = (*mpProperties)[BOND_TAU_ZERO];
    //const double& internal_friction = (*mpProperties)[BOND_INTERNAL_FRICC];
    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    double current_tangential_force_module = 0.0;

    // The [mBondedScalingFactor] is divided into two direction [0] and [1]. June, 2022
    double JointLocalElasticContactForce[2] = {0.0};

    // bond force for joint contact
    if (!failure_type) {
        //double JointSlidingLocalVel[3] = {0.0};
        //double temp_local_vel = GeometryFunctions::DotProduct(LocalRelVel, mLocalJointNormal);
        //JointSlidingLocalVel[0] =  temp_local_vel * mLocalJointNormal[0];
        //JointSlidingLocalVel[1] =  temp_local_vel * mLocalJointNormal[1];
        //LocalDeltSlidingDisp[0] = JointSlidingLocalVel[0] * r_process_info[DELTA_TIME];
        //LocalDeltSlidingDisp[1] = JointSlidingLocalVel[1] * r_process_info[DELTA_TIME];
        mAccumulatedJointTangentialLocalDisplacement[0] += LocalDeltDisp[0];
        mAccumulatedJointTangentialLocalDisplacement[1] += LocalDeltDisp[1];
        JointLocalElasticContactForce[0] -= kt_el * mAccumulatedJointTangentialLocalDisplacement[0]; // 0: first tangential
        JointLocalElasticContactForce[1] -= kt_el * mAccumulatedJointTangentialLocalDisplacement[1]; // 1: second tangential
    } else {
        mAccumulatedJointTangentialLocalDisplacement[0] += LocalRelVel[0];
        mAccumulatedJointTangentialLocalDisplacement[1] += LocalRelVel[1];
        JointLocalElasticContactForce[0] -= kt_el * mAccumulatedJointTangentialLocalDisplacement[0]; // 0: first tangential
        JointLocalElasticContactForce[1] -= kt_el * mAccumulatedJointTangentialLocalDisplacement[1]; // 1: second tangential
        current_tangential_force_module = sqrt(JointLocalElasticContactForce[0] * JointLocalElasticContactForce[0]
                                                 + JointLocalElasticContactForce[1] * JointLocalElasticContactForce[1]);
        double friction_force = 0.5 * LocalElasticContactForce[2];
        if (current_tangential_force_module > friction_force){
            if (current_tangential_force_module > 0.0){
                const double fraction = friction_force / current_tangential_force_module;
                JointLocalElasticContactForce[0] *= fraction;
                JointLocalElasticContactForce[1] *= fraction;
            }
        }
    }

    current_tangential_force_module = sqrt(JointLocalElasticContactForce[0] * JointLocalElasticContactForce[0]
                                                 + JointLocalElasticContactForce[1] * JointLocalElasticContactForce[1]);


    if (calculation_area){
        contact_tau = current_tangential_force_module / calculation_area;
    }

    LocalElasticContactForce[0] = JointLocalElasticContactForce[0];
    LocalElasticContactForce[1] = JointLocalElasticContactForce[1];

    KRATOS_CATCH("")
}

//*************************************
// Moment calculation
//*************************************

void DEM_smooth_joint::CalculateMoments(SphericContinuumParticle* element, 
                    SphericContinuumParticle* neighbor, 
                    double equiv_young, 
                    double distance, 
                    double calculation_area,
                    double LocalCoordSystem[3][3], 
                    double ElasticLocalRotationalMoment[3], 
                    double ViscoLocalRotationalMoment[3], 
                    double equiv_poisson, 
                    double indentation, 
                    double LocalElasticContactForce[3],
                    double normalLocalContactForce,
                    double GlobalElasticContactForces[3],
                    double LocalCoordSystem_2[3],
                    const int i_neighbor_count) 
{
    KRATOS_TRY

    //here we do not calculate the moment
    //according to Luc Scholtes (2012) Modelling progressive failure in fractured rock masses using a 3D discrete element method

    KRATOS_CATCH("")
}

//*************************************
// Joint failure checking
//*************************************

void DEM_smooth_joint::CheckFailure(const int i_neighbour_count, 
                                        SphericContinuumParticle* element1, 
                                        SphericContinuumParticle* element2,
                                        double& contact_sigma,
                                        double& contact_tau, 
                                        double LocalElasticContactForce[3],
                                        double ViscoDampingLocalContactForce[3],
                                        double ElasticLocalRotationalMoment[3],
                                        double ViscoLocalRotationalMoment[3]){
    
    KRATOS_TRY

    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

    if (failure_type == 0) {

        //parameters
        const double& bond_sigma_max = (*mpProperties)[BOND_SIGMA_MAX];
        //const double& bond_sigma_max_deviation = (*mpProperties)[BOND_SIGMA_MAX_DEVIATION];
        const double& bond_tau_zero = (*mpProperties)[BOND_TAU_ZERO];
        //const double& bond_tau_zero_deviation = (*mpProperties)[BOND_TAU_ZERO_DEVIATION];
        const double& bond_interanl_friction = (*mpProperties)[BOND_INTERNAL_FRICC];

        double bond_current_tau_max = bond_tau_zero;

        if (contact_sigma >= 0) {
            bond_current_tau_max += tan(bond_interanl_friction * Globals::Pi / 180.0) * contact_sigma;
        }

        if (contact_sigma < 0.0  /*break only in tension*/
                && (-1 * contact_sigma > bond_sigma_max) 
                && !(*mpProperties)[IS_UNBREAKABLE]) 
        { //for normal
            failure_type = 4; // failure in tension
            contact_sigma = 0.0;
            contact_tau = 0.0;
            LocalElasticContactForce[0] = 0.0;      
            LocalElasticContactForce[1] = 0.0;      
            LocalElasticContactForce[2] = 0.0;
        }
        else if(( std::abs(contact_tau) > bond_current_tau_max) && !(*mpProperties)[IS_UNBREAKABLE]) 
        { //for tangential 
            failure_type = 2; // failure in shear
            contact_sigma = 0.0;
            contact_tau = 0.0;
            LocalElasticContactForce[2] = 0.0;
            double current_tangential_force_module = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            double friction_force = (*mpProperties)[JOINT_FRICTION_COEFF] * LocalElasticContactForce[2];
            if (current_tangential_force_module > friction_force){
                if (current_tangential_force_module > 0.0){
                    const double fraction = friction_force / current_tangential_force_module;
                    LocalElasticContactForce[0] *= fraction;
                    LocalElasticContactForce[1] *= fraction;
                }
            }
        } 
        
    }

    KRATOS_CATCH("")    
}//CheckFailure
    
} //namespace Kratos