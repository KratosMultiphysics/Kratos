//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

// System includes
#include <algorithm>

// External includes 

// Project includes
#include "custom_utilities/circular_particle_hertzian.h"

namespace Kratos{

void CircularHertzianParticle::ComputeForcesOnCenterNode(double dt, array_1d<double, 3 >& gravity){
    double mass = mMass;
//    double n_exponent = 1.5;
    double radius = mRadius;
//    double young_mod = mYoung;
    double poisson_coef = mPoisson;
    double young_star = mYoungStar;
    double zeta = mZeta;
    double dt_inv = 1.0 / dt;
    double n_of_neighbours_inv = 1.0 / mNumberOfNeighbours;
    array_1d<double,3>& force = GetForce();
    noalias(force) = ZeroVector(3);
    force = mass * gravity;
//    DistanceIteratorType idistance = SphericHertzianParticle::mDistancesToNeighbours.begin();
    for(ParticleWeakIteratorType ineighbour = CircularHertzianParticle::mNeighbours.begin();
        ineighbour != CircularHertzianParticle::mNeighbours.end(); ineighbour++){
//        KRATOS_WATCH('In loop for adding every force');
        double other_radius = ineighbour->GetRadius();
//        double other_young = ineighbour->GetYoung();
//        double other_poisson = ineighbour->GetPoisson();
        double other_zeta = ineighbour->GetZeta();
        double other_mass = ineighbour->GetMass();
        array_1d<double,3> other_to_me_dist = this->GetPosition() - ineighbour->GetPosition();
//        distance_2 = *idistance;
        double distance_2 = other_to_me_dist[0] * other_to_me_dist[0] + other_to_me_dist[1] * other_to_me_dist[1] + other_to_me_dist[2] * other_to_me_dist[2];   //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_sum = radius + other_radius;
        double aux_1 = distance_2 - radius_sum * radius_sum;
//        KRATOS_WATCH(distance_2);
//        KRATOS_WATCH(mProximity_Tol);
        if (aux_1 < 0.0){
            double distance = sqrt(distance_2);
            double indentation = radius_sum - distance;
//            KRATOS_WATCH(distance_2);
//            KRATOS_WATCH(distance);
//            KRATOS_WATCH(indentation);
//            indentation = 0.5 * (fabs(0.00001 - indentation) + 0.00001 + indentation);
//            KRATOS_WATCH(indentation_n_exp);
//            double indentation_n_exp = sqrt(indentation * indentation * indentation);
            double other_young_star = ineighbour->GetYoungStar();
//            KRATOS_WATCH(young_star);
//            KRATOS_WATCH(other_young_star);
//            double equiv_radius_double = 2.0 / ((1.0 / radius) + (1.0 / other_radius));
            double equiv_young_star = 1.0 / ((1.0 / young_star) + (1.0 / other_young_star));
            double mod_force = 0.25 * M_PI * equiv_young_star * indentation;
//            KRATOS_WATCH(equiv_young_star);
//            KRATOS_WATCH(equiv_radius_double);
//            KRATOS_WATCH(sqrt(equiv_radius_double));
//            KRATOS_WATCH(mod_force);
            array_1d<double,3> versor = other_to_me_dist / distance;
            if (zeta + other_zeta > mProximity_Tol){
                double stiffness_quot = other_young_star * other_radius / (young_star * radius);
                double aux_2 = (1 / (1 + stiffness_quot));
                double equiv_damp_coef = poisson_coef * exp(log(aux_2));
                double equiv_damp = zeta * equiv_damp_coef * indentation;
                KRATOS_WATCH(equiv_damp);
                array_1d<double,3>& vel = GetVelocity();
                array_1d<double,3>& other_vel = ineighbour->GetVelocity();
                array_1d<double,3> other_to_me_dist_dot = vel - other_vel;
                array_1d<double,3> other_to_me_vel = (other_to_me_dist_dot[0] * versor[0] + other_to_me_dist_dot[1] * versor[1] + other_to_me_dist_dot[2] * versor[2]) * versor;
                array_1d<double,3> elastic_f = mod_force * versor;
                double aux_3 = equiv_damp * dt / mass;
                double min_mass = 0.5 * (mass + other_mass - fabs(mass - other_mass));
                double other_n_of_neighbours_inv = 1 / ineighbour->GetNumberOfNeighbours();
                double min_number_of_neighbours_inv = 0.5 * (n_of_neighbours_inv + other_n_of_neighbours_inv - fabs(n_of_neighbours_inv - other_n_of_neighbours_inv));
                if (aux_3 > 1){
                    equiv_damp = min_mass * dt_inv * min_number_of_neighbours_inv;
                    }
                array_1d<double,3> viscous_f = - equiv_damp * other_to_me_vel;
//                double elastic_mod_2 = elastic_f[0] * elastic_f[0] + elastic_f[1] * elastic_f[1] + elastic_f[2] * elastic_f[2];
//                double viscous_mod_2 = viscous_f[0] * viscous_f[0] + viscous_f[1] * viscous_f[1] + viscous_f[2] * viscous_f[2];
//                if (elastic_mod_2 > viscous_mod_2){
                noalias(force) += elastic_f + viscous_f;
                }
            else{
                array_1d<double,3> elastic_f = mod_force * versor;
                noalias(force) += elastic_f;
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
    } // CircularHertzianParticle::ComputeForcesOnCenterNode
      
}  // namespace Kratos.



