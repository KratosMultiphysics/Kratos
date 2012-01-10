//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

// System includes
#include <math.h>
#include <algorithm>

// External includes 

// Project includes
#include "custom_utilities/circular_particle.h"

namespace Kratos{

void CircularParticle::ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity){
    double mass = mMass;
    double radius = mRadius;
    double stiffness = mStiffness / radius;
    double undamp_nat_freq = sqrt(stiffness / mass);
    double critic_damp_fraction = mZeta;
    array_1d<double,3>& force = GetForce();
    noalias(force) = ZeroVector(3);
    force = mass * gravity;
//    DistanceIteratorType idistance = CircularParticle::mDistancesToNeighbours.begin();
    for(ParticleWeakIteratorType ineighbour = CircularParticle::mNeighbours.begin();
        ineighbour != CircularParticle::mNeighbours.end(); ineighbour++){
//        KRATOS_WATCH('In loop for adding every force');
        double other_mass = ineighbour->GetMass();
        double other_radius = ineighbour->GetRadius();
        double other_stiffness = ineighbour->GetStiffness() / other_radius;
        double other_undamp_nat_freq = sqrt(other_stiffness / other_mass);
        double other_critic_damp_fraction = ineighbour->GetZeta() * other_stiffness;
        
        array_1d<double,3> other_to_me_dist = this->GetPosition() - ineighbour->GetPosition();
//        distance_2 = *idistance;
        double distance_2 = other_to_me_dist[0] * other_to_me_dist[0] + other_to_me_dist[1] * other_to_me_dist[1] + other_to_me_dist[2] * other_to_me_dist[2];   //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_sum = radius + other_radius;
        double aux_1 = distance_2 - radius_sum * radius_sum;
//        KRATOS_WATCH(distance_2);
//        KRATOS_WATCH(mProximity_Tol);
        if (aux_1 < mProximity_Tol){
            double distance = sqrt(distance_2);
            array_1d<double,3> versor = other_to_me_dist / distance;
            double equiv_stiffness = stiffness * other_stiffness / (other_stiffness + stiffness);
            double indentation = radius_sum - distance;
            double normal_elastic_force = equiv_stiffness * indentation;
            if (critic_damp_fraction + other_critic_damp_fraction > mProximity_Tol){
                double reduced_mass = mass * other_mass / (mass + other_mass);
                double critical_damping = 2 * sqrt(reduced_mass * equiv_stiffness);
                double damping_coef = critic_damp_fraction * critical_damping;
//                double equiv_damp = 2 / ((1 / (undamp_nat_freq * critic_damp_fraction)) + (1 / (other_undamp_nat_freq * other_critic_damp_fraction)));
                array_1d<double,3>& vel = GetVelocity();
                array_1d<double,3>& other_vel = ineighbour->GetVelocity();
                array_1d<double,3> other_to_me_dist_dot = vel - other_vel;
                double normal_other_to_me_celerity = (other_to_me_dist_dot[0] * versor[0] + other_to_me_dist_dot[1] * versor[1] + other_to_me_dist_dot[2] * versor[2]);
//                array_1d<double,3> other_to_me_vel = other_to_me_dist_dot_proj * versor;
//                double aux_2 = damping_coef * dt / mass;
//                double min_mass = 0.5 * (mass + other_mass - abs(mass - other_mass));
//                double max_frequence_2 = 0.5 * (stiffness / mass + other_stiffness / other_mass + abs(stiffness / mass - other_stiffness / other_mass));
//                double other_n_of_neighbours_inv = 1 / ineighbour->GetNumberOfNeighbours();
//                double min_number_of_neighbours_inv = 0.5 * (n_of_neighbours_inv + other_n_of_neighbours_inv - abs(n_of_neighbours_inv - other_n_of_neighbours_inv));
//                if (aux_2 > 1){
//                    equiv_damp = min_mass * dt_inv * min_number_of_neighbours_inv;
//                    }
                // Computing maximum to avoid "sticking"
                double normal_viscous_f = 0.5 * (-normal_elastic_force + damping_coef * normal_other_to_me_celerity + fabs(-normal_elastic_force - damping_coef * normal_other_to_me_celerity));
                noalias(force) += (normal_elastic_force + normal_viscous_f) * versor;
                }
            else{
                noalias(force) += normal_elastic_force * versor;
                }
//            KRATOS_WATCH(distance);
//            KRATOS_WATCH(equiv_stiffness);
//            KRATOS_WATCH(aux_2);
//            KRATOS_WATCH('Contribution force of particle:');
//            KRATOS_WATCH(ineighbour->Id());
//            KRATOS_WATCH(contribution);
            }
//      idistance++;
        }
//    KRATOS_WATCH("TOTAL FORCE on particle");
//    KRATOS_WATCH(this->Id());
//    KRATOS_WATCH(this->GetVelocity());
//    KRATOS_WATCH((*this)[0]);
//    KRATOS_WATCH((*this)[1]);
//    KRATOS_WATCH(force);
    } // CircularParticle::ComputeForcesOnCenterNode
      
}  // namespace Kratos.



