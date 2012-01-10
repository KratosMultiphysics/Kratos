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
#include "custom_utilities/spheric_particle_hertzian.h"

namespace Kratos{

void SphericHertzianParticle::ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity){
    double mass = mMass;
    double n_exponent = 1.5;
    double radius = mRadius;
    double young_star = mYoungStar;
    double critic_damp_fraction = mZeta;
    array_1d<double,3>& force = GetForce();
    noalias(force) = ZeroVector(3);
    force = mass * gravity;
//    DistanceIteratorType idistance = SphericHertzianParticle::mDistancesToNeighbours.begin();
    for(ParticleWeakIteratorType ineighbour = SphericHertzianParticle::mNeighbours.begin();
        ineighbour != SphericHertzianParticle::mNeighbours.end(); ineighbour++){
//        KRATOS_WATCH('In loop for adding every force');
        double other_radius = ineighbour->GetRadius();
        double other_critic_damp_fraction = ineighbour->GetZeta();
        double other_mass = ineighbour->GetMass();
        array_1d<double,3> other_to_me_vect = this->GetPosition() - ineighbour->GetPosition();
        double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
        double radius_sum = radius + other_radius;
        double aux_1 = distance_2 - radius_sum * radius_sum;

        if (aux_1 < 0.0){
            double distance = sqrt(distance_2);
            double indentation = radius_sum - distance;
            //double indentation_n_exp = exp(n_exponent * log(indentation));
            //double indentation_exp_quarter = exp(0.25 * log(indentation));
            double indentation_n_exp = pow(indentation, n_exponent);
            double indentation_exp_quarter = pow(indentation, 0.16666666666666666666667 * n_exponent);
            double other_young_star = ineighbour->GetYoungStar();
            double equiv_radius = radius * other_radius / (radius + other_radius);
            double equiv_young_star = young_star * other_young_star / (young_star + other_young_star);
            double eff_stiffness = 1.333333333333333333333333 * equiv_young_star * sqrt(equiv_radius);
            double normal_elastic_force = eff_stiffness * indentation_n_exp;
            array_1d<double,3> other_to_me_versor = other_to_me_vect / distance;
            if (critic_damp_fraction + other_critic_damp_fraction > mProximity_Tol){
                double reduced_mass = mass * other_mass / (mass + other_mass);
                double critical_damping_by_half = sqrt(reduced_mass * eff_stiffness);
                double damping_coef = critic_damp_fraction * critical_damping_by_half;
//                double stiffness_quot = other_young_star * other_radius / (young_star * radius);
//                double aux_2 = (1 / (1 + stiffness_quot));
//                double equiv_damp_coef = poisson_coef * exp((n_exponent + 1) * log(aux_2));
//                double equiv_damp = critic_damp_fraction * equiv_damp_coef * indentation_n_exp;
                array_1d<double,3>& vel = GetVelocity();
                array_1d<double,3>& other_vel = ineighbour->GetVelocity();
                array_1d<double,3> other_to_me_vect_dot = vel - other_vel;
                double separation_celerity = (other_to_me_vect_dot[0] * other_to_me_versor[0] + other_to_me_vect_dot[1] * other_to_me_versor[1] + other_to_me_vect_dot[2] * other_to_me_versor[2]);
                double normal_viscous_f = -damping_coef * indentation_exp_quarter * separation_celerity;
#if STICKING_ALLOWED == 0
                if(separation_celerity > 0.0){
                    double normal_viscous_f = 0.5 * (-normal_elastic_force + normal_viscous_f - fabs(-normal_elastic_force - normal_viscous_f));
                    }
#endif
                noalias(force) += (normal_elastic_force + normal_viscous_f) * other_to_me_versor;
                }
            else{
                noalias(force) += normal_elastic_force * other_to_me_versor;
                }
//            KRATOS_WATCH(this->Id());
//            KRATOS_WATCH(ineighbour->Id());
//            KRATOS_WATCH(distance);
//            KRATOS_WATCH(equiv_stiffness);
//            KRATOS_WATCH(aux_2);
//            KRATOS_WATCH('Contribution force of particle:');

//            KRATOS_WATCH(force);
            }
//      idistance++;
        }
//    KRATOS_WATCH("TOTAL FORCE on particle");
//    KRATOS_WATCH(this->Id());
//    KRATOS_WATCH(this->GetVelocity());
//    KRATOS_WATCH((*this)[0]);
//    KRATOS_WATCH((*this)[1]);
//    KRATOS_WATCH(force);
    } // SphericHertzianParticle::ComputeForcesOnCenterNode
      
}  // namespace Kratos.



