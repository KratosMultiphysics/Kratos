//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//
#define STICKING_ALLOWED 0
// System includes
#include <math.h>
#include <algorithm>

// External includes 

// Project includes
#include "custom_utilities/spheric_particle.h"

namespace Kratos{

void SphericParticle::ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity){

    double mass = mMass;
    double radius = mRadius;
    double stiffness = mStiffness / radius;
    double critic_damp_fraction = mZeta;
    array_1d<double,3>& force = GetForce();
    noalias(force) = ZeroVector(3);
    force = mass * gravity;
    for(ParticleWeakIteratorType ineighbour = SphericParticle::mNeighbours.begin();
        ineighbour != SphericParticle::mNeighbours.end(); ineighbour++){
        double other_mass = ineighbour->GetMass();
        double other_radius = ineighbour->GetRadius();
        double other_stiffness = ineighbour->GetStiffness() / other_radius;
        double other_critic_damp_fraction = ineighbour->GetZeta();       
        array_1d<double,3> other_to_me_vect = this->GetPosition() - ineighbour->GetPosition();
        double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
        double radius_sum = radius + other_radius;
        double aux_1 = distance_2 - radius_sum * radius_sum;

        if (aux_1 < mProximity_Tol){
            double distance = sqrt(distance_2);
            double equiv_stiffness = stiffness * other_stiffness / (other_stiffness + stiffness);
            double indentation = radius_sum - distance;
            double normal_elastic_force = equiv_stiffness * indentation;
            array_1d<double,3> other_to_me_versor = other_to_me_vect / distance;
            if (critic_damp_fraction + other_critic_damp_fraction > mProximity_Tol){
                double reduced_mass = mass * other_mass / (mass + other_mass);
                double critical_damping = 2.0 * sqrt(reduced_mass * equiv_stiffness);
                double damping_coef = critic_damp_fraction * critical_damping;
                array_1d<double,3>& vel = GetVelocity();
                array_1d<double,3>& other_vel = ineighbour->GetVelocity();
                array_1d<double,3> other_to_me_vect_dot = vel - other_vel;
                double separation_celerity = other_to_me_vect_dot[0] * other_to_me_versor[0] + other_to_me_vect_dot[1] * other_to_me_versor[0] + other_to_me_vect_dot[2] * other_to_me_versor[0];
                // Computing maximum to avoid "sticking"
                double normal_viscous_f = -damping_coef * separation_celerity;
#if STICKING_ALLOWED == 0
                if(separation_celerity > 0.0){
                    double normal_viscous_f = 0.5 * (-normal_elastic_force + normal_viscous_f + fabs(-normal_elastic_force - normal_viscous_f)); // The maximum between the minus elastic force and the vioscous force
                    }
#endif
                noalias(force) += (normal_elastic_force + normal_viscous_f) * other_to_me_versor;
//                KRATOS_WATCH(distance);
//                KRATOS_WATCH(separation_celerity);
//                KRATOS_WATCH(normal_elastic_force);
//                KRATOS_WATCH(normal_viscous_f);
                }
            else{
                noalias(force) += normal_elastic_force * other_to_me_versor;
                }               
    //            KRATOS_WATCH(aux_2);
    //            KRATOS_WATCH('Contribution force of particle:');
    //            KRATOS_WATCH(ineighbour->Id());
    //            KRATOS_WATCH(contribution);
            }
    //      idistance++;
        }

    } // SphericParticle::ComputeForcesOnCenterNode
      
}  // namespace Kratos.


