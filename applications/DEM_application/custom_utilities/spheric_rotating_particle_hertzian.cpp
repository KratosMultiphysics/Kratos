//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//
#define STICKING_ALLOWED 0
// System includes
#include <algorithm>

// External includes

// Project includes
#include "custom_utilities/spheric_rotating_particle_hertzian.h"

namespace Kratos
{

void SphericRotatingHertzianParticle::ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity){

    double n_exponent = 1.5;
    double mass = mMass;
    double radius = mRadius;
    double young_star = mYoungStar;
    double stiffness = mNormalStiffness / radius;
    double tang_stiffness = mTangentialStiffness;
    double critic_damp_fraction = mZeta;
    double static_friction_coef = mStaticFriction;
    double dynamic_friction_coef = mDynamicFriction;
    array_1d<double, 3 > & force = GetForce();
    array_1d<double, 3 > & moment = GetMoment();
    noalias(force) = ZeroVector(3);
    noalias(moment) = ZeroVector(3);
    force = mass * gravity;
    ForceNormsVectorType& max_tang_forces = mMaxTangForces;
    max_tang_forces.clear();
    BoolSlidingContactsVectorType& sliding_occurrences = mSlidingFlagsOfNeighbours;
    sliding_occurrences.clear();
    TangDisplacementsIteratorType idisplacement = mTangentialDisplacementsOfNeighbours.begin();
    for (ParticleWeakIteratorType ineighbour = SphericRotatingHertzianParticle::mContactingNeighbours.begin();
            ineighbour != SphericRotatingHertzianParticle::mContactingNeighbours.end(); ineighbour++){
        array_1d<double, 3 > tang_displacement = *(idisplacement.base());
        double other_mass = ineighbour->GetMass();
        double other_radius = ineighbour->GetRadius();
        double other_young_star = ineighbour->GetYoungStar();
        double other_tang_stiffness = ineighbour->GetTangentialStiffness() / other_radius;
        double other_critic_damp_fraction = ineighbour->GetZeta();
        array_1d<double, 3 > other_to_me_vect = this->GetPosition() - ineighbour->GetPosition();
        array_1d<double, 3 > other_to_me_vect_incr = other_to_me_vect - this->GetOldPosition() + ineighbour->GetOldPosition();
        double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
        double radius_sum = radius + other_radius;
        double equiv_radius = radius * other_radius / (radius + other_radius);
        double distance = sqrt(distance_2);
        array_1d<double, 3 > other_to_me_versor = other_to_me_vect / distance;
        double equiv_young_star = young_star * other_young_star / (young_star + other_young_star);
        double eff_stiffness = 1.333333333333333333333333 * equiv_young_star * sqrt(equiv_radius);
        double equiv_tang_stiffness = tang_stiffness * other_tang_stiffness / (tang_stiffness + other_tang_stiffness);
        double indentation = radius_sum - distance;
        double indentation_n_exp = pow(indentation, n_exponent);
        double indentation_exp_quarter = pow(indentation, 0.16666666666666666666667 * n_exponent);
        double normal_force_norm = eff_stiffness * indentation_n_exp;//Minus the viscous part (Elastic case)
        array_1d<double, 3 > tang_force = -equiv_tang_stiffness * tang_displacement;//Minus the viscous part (Elastic case)
        double calculated_tang_force_norm = sqrt(tang_force[0] * tang_force[0] + tang_force[1] * tang_force[1] + tang_force[2] * tang_force[2]);//Elastic case
        if (critic_damp_fraction + other_critic_damp_fraction > mProximity_Tol){
            double reduced_mass = mass * other_mass / (mass + other_mass);
            double critical_damping = 2.0 * sqrt(reduced_mass * eff_stiffness);
            double normal_damping_coef = critic_damp_fraction * critical_damping;
            double tang_damping_coef = normal_damping_coef; //Tsuji et al. (1992)
            array_1d<double, 3 > & vel = GetVelocity();
            array_1d<double, 3 > & ang_vel = GetAngularVelocity();
            array_1d<double, 3 > & other_ang_vel = ineighbour->GetAngularVelocity();
            array_1d<double, 3 > & other_vel = ineighbour->GetVelocity();
            array_1d<double, 3 > other_to_me_vect_dot = vel - other_vel;
            double separation_celerity = other_to_me_vect_dot[0] * other_to_me_versor[0] + other_to_me_vect_dot[1] * other_to_me_versor[0] + other_to_me_vect_dot[2] * other_to_me_versor[0];
            array_1d<double, 3 > tang_velocity_of_me_from_other = other_to_me_vect_dot - separation_celerity * other_to_me_versor;
            double tang_celerity_of_me_from_other = sqrt(tang_velocity_of_me_from_other[0] * tang_velocity_of_me_from_other[0] + tang_velocity_of_me_from_other[1] * tang_velocity_of_me_from_other[1] + tang_velocity_of_me_from_other[2] * tang_velocity_of_me_from_other[2]); ;
            // Computing maximum to avoid "sticking"
            double normal_viscous_f = -normal_damping_coef * separation_celerity;
            double tang_viscous_f = -tang_damping_coef * tang_celerity_of_me_from_other;
#if STICKING_ALLOWED == 0
            if (separation_celerity > 0.0){
                double normal_viscous_f = 0.5 * (-normal_force_norm + normal_viscous_f + fabs(-normal_force_norm - normal_viscous_f)); // The maximum between the minus elastic force and the vioscous force
            }
#endif
            normal_force_norm = normal_force_norm + normal_viscous_f;
//#if STICKING_ALLOWED == 0
//            normal_force_norm = fabs(normal_force_norm + normal_viscous_f);
//#endif
            array_1d<double, 3 > tang_versor = tang_displacement;
            array_1d<double, 3 > contribution_of_spin_to_vel, other_contribution_of_spin_to_vel;
            contribution_of_spin_to_vel[0] = radius * (ang_vel[2] * other_to_me_versor[1] - ang_vel[1] * other_to_me_versor[2]);
            contribution_of_spin_to_vel[1] = radius * (ang_vel[0] * other_to_me_versor[2] - ang_vel[2] * other_to_me_versor[0]);
            contribution_of_spin_to_vel[2] = radius * (ang_vel[1] * other_to_me_versor[0] - ang_vel[0] * other_to_me_versor[1]);
            other_contribution_of_spin_to_vel[0] = other_radius * (other_ang_vel[1] * other_to_me_versor[2] - other_ang_vel[2] * other_to_me_versor[1]);
            other_contribution_of_spin_to_vel[1] = other_radius * (other_ang_vel[2] * other_to_me_versor[0] - other_ang_vel[0] * other_to_me_versor[2]);
            other_contribution_of_spin_to_vel[2] = other_radius * (other_ang_vel[0] * other_to_me_versor[1] - other_ang_vel[1] * other_to_me_versor[0]);
            array_1d<double, 3 > difference_contact_velocities = other_to_me_vect_dot - contribution_of_spin_to_vel - other_contribution_of_spin_to_vel;
            array_1d<double, 3 > tang_vel_vector = difference_contact_velocities - (difference_contact_velocities[0] * other_to_me_versor[0] + difference_contact_velocities[1] * other_to_me_versor[1] + difference_contact_velocities[2] * other_to_me_versor[2]) * other_to_me_versor;
            double tang_vel_norm = sqrt(tang_vel_vector[0] * tang_vel_vector[0] + tang_vel_vector[1] * tang_vel_vector[1] + tang_vel_vector[2] * tang_vel_vector[2]);
            if(tang_vel_norm > mProximity_Tol){
                array_1d<double, 3 > tang_vel_versor = tang_vel_vector / tang_vel_norm;
                double tang_displ_norm_2 = tang_displacement[0] * tang_displacement[0] + tang_displacement[1] * tang_displacement[1] + tang_displacement[2] * tang_displacement[2];
                double tang_vel_versor_proj_on_displ_direction;
                if(tang_displ_norm_2 > mProximity_Tol)
                    tang_vel_versor_proj_on_displ_direction = (tang_vel_vector[0] * tang_displacement[0] + tang_vel_vector[1] * tang_displacement[1] + tang_vel_vector[2] * tang_displacement[2]) / sqrt(tang_displ_norm_2);
                else{
                    tang_vel_versor_proj_on_displ_direction = 0.0;
                }
                double tang_viscous_force = tang_damping_coef * tang_vel_versor_proj_on_displ_direction * tang_vel_norm;
                tang_force = tang_force + tang_viscous_force * tang_vel_versor;
                calculated_tang_force_norm = sqrt(tang_force[0] * tang_force[0] + tang_force[1] * tang_force[1] + tang_force[2] * tang_force[2]);
                double aux_positive_in_loading = tang_vel_versor[0] * tang_displacement[0] + tang_vel_versor[1] * tang_displacement[0] + tang_vel_versor[2] * tang_displacement[0];
#if STICKING_ALLOWED == 0
                if(aux_positive_in_loading < 0.0){
                    double tang_viscous_force = 0.5 * (-calculated_tang_force_norm + tang_viscous_force + fabs(-calculated_tang_force_norm - tang_viscous_force)); // The maximum between the minus elastic force and the vioscous force
                }
            }
#endif
        }
        double max_tang_force_norm = static_friction_coef * normal_force_norm;
        max_tang_forces.push_back(max_tang_force_norm); //Needed for SphericRotatingHertzianParticle::ComputeNewTangentialDisplacements()
        double tang_force_norm = 0.5 * (calculated_tang_force_norm + max_tang_force_norm - fabs(calculated_tang_force_norm - max_tang_force_norm));
        if (calculated_tang_force_norm > max_tang_force_norm){
            sliding_occurrences.push_back(true);
        }
        else{
            sliding_occurrences.push_back(false);
        }
        tang_force = tang_force * tang_force_norm / (calculated_tang_force_norm + 0.01);
        noalias(force) += normal_force_norm * other_to_me_versor + tang_force;
        moment[0] += radius * (tang_force[1] * other_to_me_versor[2] - tang_force[2] * other_to_me_versor[1]);
        moment[1] += radius * (tang_force[2] * other_to_me_versor[0] - tang_force[0] * other_to_me_versor[2]);
        moment[2] += radius * (tang_force[0] * other_to_me_versor[1] - tang_force[1] * other_to_me_versor[0]);
        ++idisplacement;
    }
} // SphericRotatingHertzianParticle::ComputeForcesOnCenterNode

void SphericRotatingHertzianParticle::ComputeNewTangentialDisplacements(){
    int n_neigbours = SphericRotatingHertzianParticle::GetNumberOfNeighbours();
    int n_contacts = SphericRotatingHertzianParticle::GetNumberOfContactingNeighbours();
    if (n_neigbours > 0 || (n_neigbours < 1 && n_contacts > 0)){
        double radius = mRadius;
        array_1d<double, 3 > rotation_incr = this->GetRotation() - this->GetOldRotation();
        double phi = rotation_incr[0];
        double theta = rotation_incr[1];
        double omega = rotation_incr[2];
        double cp = cos(phi);
        double ct = cos(theta);
        double co = cos(omega);
        double sp = sin(phi);
        double st = sin(theta);
        double so = sin(omega);
        double rot_mat_11 = ct * co;
        double rot_mat_12 = sp * st * co - cp * so;
        double rot_mat_13 = sp * so + cp * st * co;
        double rot_mat_21 = ct * so;
        double rot_mat_22 = sp * st * so + cp * co;
        double rot_mat_23 = cp * st * so - sp * co;
        double rot_mat_31 = -st;
        double rot_mat_32 = sp * ct;
        double rot_mat_33 = cp * ct;
        ForceNormsIteratorType inormalf = SphericRotatingHertzianParticle::mMaxTangForces.begin();
        BoolSlidingContactsIteratorType islide = SphericRotatingHertzianParticle::mSlidingFlagsOfNeighbours.begin();
        TangDisplacementsIteratorType idisplacement = SphericRotatingHertzianParticle::mTangentialDisplacementsOfNeighbours.begin();
        for (ParticleWeakIteratorType ineighbour = SphericRotatingHertzianParticle::mContactingNeighbours.begin();
        ineighbour != SphericRotatingHertzianParticle::mContactingNeighbours.end(); ineighbour++){
            //array_1d<double, 3 > tang_displacement = *(idisplacement.base());
            double other_radius = ineighbour->GetRadius();
            array_1d<double, 3 >& tang_displ = *idisplacement;
            array_1d<double, 3 > other_to_me_vect = this->GetPosition() - ineighbour->GetPosition();
            array_1d<double, 3 > other_to_me_vect_incr = other_to_me_vect - this->GetOldPosition() + ineighbour->GetOldPosition();
            array_1d<double, 3 > other_rotation_incr = ineighbour->GetRotation() - ineighbour->GetOldRotation();
            double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
            double distance = sqrt(distance_2);
            array_1d<double, 3 > other_to_me_versor = other_to_me_vect / distance;
            array_1d<double, 3 > contribution_of_ang_incr_to_tang_displ, other_contribution_of_ang_incr_to_tang_displ;
            contribution_of_ang_incr_to_tang_displ[0] = radius * (rotation_incr[2] * other_to_me_versor[1] - rotation_incr[1] * other_to_me_versor[2]);
            contribution_of_ang_incr_to_tang_displ[1] = radius * (rotation_incr[0] * other_to_me_versor[2] - rotation_incr[2] * other_to_me_versor[0]);
            contribution_of_ang_incr_to_tang_displ[2] = radius * (rotation_incr[1] * other_to_me_versor[0] - rotation_incr[0] * other_to_me_versor[1]);
            other_contribution_of_ang_incr_to_tang_displ[0] = other_radius * (other_rotation_incr[1] * other_to_me_versor[2] - other_rotation_incr[2] * other_to_me_versor[1]);
            other_contribution_of_ang_incr_to_tang_displ[1] = other_radius * (other_rotation_incr[2] * other_to_me_versor[0] - other_rotation_incr[0] * other_to_me_versor[2]);
            other_contribution_of_ang_incr_to_tang_displ[2] = other_radius * (other_rotation_incr[0] * other_to_me_versor[1] - other_rotation_incr[1] * other_to_me_versor[0]);
            array_1d<double, 3 > difference_contact_displ = other_to_me_vect_incr - contribution_of_ang_incr_to_tang_displ - other_contribution_of_ang_incr_to_tang_displ;
            array_1d<double, 3 > tang_displacement_incr = difference_contact_displ - (difference_contact_displ[0] * other_to_me_versor[0] + difference_contact_displ[1] * other_to_me_versor[1] + difference_contact_displ[2] * other_to_me_versor[2]) * other_to_me_versor;
            array_1d<double, 3 > calculated_tang_displ = tang_displ + tang_displacement_incr;
            double old_tang_displ_norm_2 = tang_displ[0] * tang_displ[0] + tang_displ[1] * tang_displ[1] + tang_displ[2] * tang_displ[2];
            double calculated_tang_displ_norm_2 = calculated_tang_displ[0] * calculated_tang_displ[0] + calculated_tang_displ[1] * calculated_tang_displ[1] + calculated_tang_displ[2] * calculated_tang_displ[2];
            double max_tang_force_norm = *inormalf;
            double max_tang_displ_2 = max_tang_force_norm * max_tang_force_norm / (mTangentialStiffness * mTangentialStiffness);
            noalias(tang_displ) += tang_displacement_incr;
            if(max_tang_displ_2 < calculated_tang_displ_norm_2 && *islide){
                noalias(tang_displ) =  tang_displ * sqrt(max_tang_displ_2 / calculated_tang_displ_norm_2);
            }
            double td1 = tang_displ[0];
            double td2 = tang_displ[1];
            double td3 = tang_displ[2];
            tang_displ[0] = rot_mat_11 * td1 + rot_mat_12 * td2 + rot_mat_13 * td3;
            tang_displ[1] = rot_mat_21 * td1 + rot_mat_22 * td2 + rot_mat_23 * td3;
            tang_displ[2] = rot_mat_31 * td1 + rot_mat_32 * td2 + rot_mat_33 * td3;
            ++idisplacement;
            ++islide;
            ++inormalf;
        }
    }
} // SphericRotatingHertzianParticle::ComputeNewTangentialDisplacements

void SphericRotatingHertzianParticle::UpdateContactsList(){
    int n_neigbours = SphericRotatingHertzianParticle::GetNumberOfNeighbours();
    int n_contacts = SphericRotatingHertzianParticle::GetNumberOfContactingNeighbours();
    if (n_neigbours > 0 || (n_neigbours < 1 && n_contacts > 0)){
        array_1d<double, 3 > init_tang_displ = ZeroVector(3);
        double radius = mRadius;
        ParticleWeakVectorType& contact_neighbours = SphericRotatingHertzianParticle::mContactingNeighbours;
        ParticleWeakVectorType old_contact_neighbours;
        TangDisplacementsVectorType& tangential_displacements = SphericRotatingHertzianParticle::mTangentialDisplacementsOfNeighbours;
        TangDisplacementsVectorType old_tangential_displacements;
        old_contact_neighbours.reserve(SphericRotatingHertzianParticle::mNeighbours.size());
        old_tangential_displacements.reserve(SphericRotatingHertzianParticle::mNeighbours.size());
        old_contact_neighbours.swap(contact_neighbours);
        old_tangential_displacements.swap(tangential_displacements);
        for (ParticleWeakIteratorType ineighbour = mNeighbours.begin();
        ineighbour != mNeighbours.end(); ++ineighbour) {
            array_1d<double, 3 > other_to_me_vect = this->GetPosition() - ineighbour->GetPosition();
            double other_radius = ineighbour->GetRadius();
            double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
            double radius_sum = radius + other_radius;
            double aux = distance_2 - radius_sum * radius_sum;
            if (aux < 0.0){
                contact_neighbours.push_back(*(ineighbour.base()));
                TangDisplacementsIteratorType  idisplacement = old_tangential_displacements.begin();
                bool found = false;
                for (ParticleWeakIteratorType icontactingneighbour = old_contact_neighbours.begin();
                icontactingneighbour != old_contact_neighbours.end(); icontactingneighbour++){
                    if(icontactingneighbour->GetPointerToCenterNode() == ineighbour->GetPointerToCenterNode()){
                        found = true;
                        tangential_displacements.push_back(*(idisplacement.base()));
                        break;
                    }
                    idisplacement++;
                }
                if(found == false){
                    tangential_displacements.push_back(init_tang_displ);
                }
            }
        }
    }
}//SphericRotatingHertzianParticle::UpdateContactsList

} // namespace Kratos.


