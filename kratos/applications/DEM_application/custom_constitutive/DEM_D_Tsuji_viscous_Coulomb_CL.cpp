// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_D_Tsuji_viscous_Coulomb_CL.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Tsuji_viscous_Coulomb::Initialize(const ProcessInfo& rCurrentProcessInfo){
        
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Tsuji_viscous_Coulomb::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Tsuji_viscous_Coulomb(*this));
        return p_clone;
    }

    void DEM_D_Tsuji_viscous_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << " Assigning DEM_D_Tsuji_viscous_Coulomb to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }                    
        
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////////NEW IMPLEMENTATION/////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateGamma( const double e, double & gamma){
    
        if(e==0.0) {
            gamma = 20.0;//Approximate value for restit. coeff of 0.001
            return;
        }
        const double h1  = -6.918798;
        const double h2  = -16.41105;
        const double h3  =  146.8049;
        const double h4  = -796.4559;
        const double h5  =  2928.711;
        const double h6  = -7206.864;
        const double h7  =  11494.29;
        const double h8  = -11342.18;
        const double h9  =  6276.757;
        const double h10 = -1489.915;        
        
        const double alpha = e*(h1+e*(h2+e*(h3+e*(h4+e*(h5+e*(h6+e*(h7+e*(h8+e*(h9+e*h10)))))))));
        
        gamma = sqrt(1.0/(1.0 - (1.0+e)*(1.0+e) * exp(alpha)) - 1.0);     
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::InitializeContact(SphericParticle* const element1,
                                                 SphericParticle* const element2) {
        //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;
        
        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear_modulus = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);                   
        
        //Normal and Tangent elastic constants
        mKn = 1.3333333333333333 * equiv_young * sqrt(equiv_radius);
        mKt = 8.0 * equiv_shear_modulus * sqrt(equiv_radius);
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateForces(const double OldLocalContactForce[3],
                                                             double LocalElasticContactForce[3],
                                                             double LocalDeltDisp[3],
                                                             double LocalRelVel[3],            
                                                             double indentation,
                                                             double previous_indentation,
                                                             double ViscoDampingLocalContactForce[3],
                                                             double& cohesive_force,
                                                             SphericParticle* element1,
                                                             SphericParticle* element2) {
        
        bool sliding = false;
        
        InitializeContact(element1, element2);
        
        LocalElasticContactForce[2]  = CalculateNormalForce(indentation, element1, element2);
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2);                                                      
        
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, sliding, element1, element2, indentation);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForce(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element1, element2, indentation, previous_indentation);
        
        if (sliding) ViscoDampingLocalContactForce[0] = ViscoDampingLocalContactForce[1] = 0.0;                
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2,
                                                                double indentation) {                                        
        
        //Get equivalent mass
        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
                
        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);        
        
        //Get equivalent ln of restitution coefficient
        const double my_ln_of_restit_coeff    = element1->GetLnOfRestitCoeff();
        const double other_ln_of_restit_coeff = element2->GetLnOfRestitCoeff();
        
        double my_coefficient_of_restitution;
        if (my_ln_of_restit_coeff > 0.0) my_coefficient_of_restitution = 0.0;
        else my_coefficient_of_restitution = exp(my_ln_of_restit_coeff);

        double other_coefficient_of_restitution;
        if (other_ln_of_restit_coeff > 0.0) other_coefficient_of_restitution = 0.0;        
        else other_coefficient_of_restitution = exp(other_ln_of_restit_coeff);

        const double equiv_coefficient_of_restitution = 0.5 * (my_coefficient_of_restitution + other_coefficient_of_restitution);        
        
        double normal_damping_coefficient;
        
        normal_damping_coefficient = (0.833333333333 * equiv_coefficient_of_restitution * equiv_coefficient_of_restitution - 2.25 * equiv_coefficient_of_restitution + 1.41666666666666) * sqrt(equiv_mass * mKn) * sqrt(sqrt(indentation));
        
        const double tangential_damping_coefficient = normal_damping_coefficient;
        
        double gamma = 0.0;
        CalculateGamma(equiv_coefficient_of_restitution, gamma);        
        normal_damping_coefficient = 2.0 * gamma * sqrt(equiv_mass * mKn) * sqrt(sqrt(indentation));
        
        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient     * LocalRelVel[2];
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateTangentialForce(const double normal_contact_force,
                                                                const double OldLocalContactForce[3],
                                                                double LocalElasticContactForce[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                const double LocalDeltDisp[3],            
                                                                bool& sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2,
                                                                double indentation,
                                                                double previous_indentation) { 
        
        LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * sqrt(indentation) * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * sqrt(indentation) * LocalDeltDisp[1];
        
        const double my_tg_of_friction_angle    = element1->GetTgOfFrictionAngle();
        const double other_tg_of_friction_angle = element2->GetTgOfFrictionAngle();
        const double equiv_tg_of_fri_ang        = 0.5 * (my_tg_of_friction_angle + other_tg_of_friction_angle);                
        
        const double MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_fri_ang;
        
        double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];
        
        double ActualTotalShearForce   = sqrt(tangential_contact_force_0  * tangential_contact_force_0  + tangential_contact_force_1  * tangential_contact_force_1);
                
        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {
            
            const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            
            if (MaximumAdmisibleShearForce < ActualElasticShearForce) {
                LocalElasticContactForce[0] = LocalElasticContactForce[0] * MaximumAdmisibleShearForce / ActualElasticShearForce;
                LocalElasticContactForce[1] = LocalElasticContactForce[1] * MaximumAdmisibleShearForce / ActualElasticShearForce;
                ViscoDampingLocalContactForce[0] = 0.0;
                ViscoDampingLocalContactForce[1] = 0.0;
            }
            
            else {
                
                const double ActualViscousShearForce = MaximumAdmisibleShearForce - ActualElasticShearForce;
                double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0]\
                                                        + ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);
        
                ViscoDampingLocalContactForce[0] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                ViscoDampingLocalContactForce[1] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
            }
                    
            sliding = true;   
        } 
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta) {
        
        const double my_radius  = element->GetRadius(); //Get equivalent Radius
        const double my_young   = element->GetYoung(); //Get equivalent Young's Modulus
        const double my_poisson = element->GetPoisson();
                
        //Normal and Tangent elastic constants
        mKn = (1.33333333 * my_young / (1.0 - my_poisson * my_poisson)) * sqrt(my_radius - ini_delta); //TODO: CHECK THIS 
        mKt = 4.0 * sqrt(my_radius) * my_young / ((1.0 + my_poisson) * (2.0 - my_poisson)); //LATER ON THIS VALUE WILL BE MULTIPLIED BY SQRT(INDENTATION)                       
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateForcesWithFEM(const double OldLocalContactForce[3],
                                                             double LocalElasticContactForce[3],
                                                             double LocalDeltDisp[3],
                                                             double LocalRelVel[3],            
                                                             double indentation,
                                                             double previous_indentation,
                                                             double ViscoDampingLocalContactForce[3],
                                                             double& cohesive_force,
                                                             SphericParticle* const element,
                                                             DEMWall* const wall) {
        bool sliding = false;
        
        InitializeContactWithFEM(element, wall);
        
        LocalElasticContactForce[2]  = CalculateNormalForceWithFEM(indentation, element, wall);
        cohesive_force               = CalculateCohesiveNormalForceWithFEM(element, wall);                                                      
        
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, sliding, element, wall, indentation);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForceWithFEM(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element, wall, indentation, previous_indentation);
        
        if (sliding) ViscoDampingLocalContactForce[0] = ViscoDampingLocalContactForce[1] = 0.0;
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateTangentialForceWithFEM(const double normal_contact_force,
                                                                      const double OldLocalContactForce[3],
                                                                      double LocalElasticContactForce[3],
                                                                      double ViscoDampingLocalContactForce[3],
                                                                      const double LocalDeltDisp[3],            
                                                                      bool& sliding,
                                                                      SphericParticle* const element,
                                                                      DEMWall* const wall,
                                                                      double indentation,
                                                                      double previous_indentation) {
                                                                        
        LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * sqrt(indentation) * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * sqrt(indentation) * LocalDeltDisp[1];
        
        const double WallBallFriction = wall->mTgOfFrictionAngle;
        
        const double MaximumAdmisibleShearForce = normal_contact_force * WallBallFriction;
        
        const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];
        
        const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);
                
        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {
            
            const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            
            if (MaximumAdmisibleShearForce < ActualElasticShearForce) {
                LocalElasticContactForce[0] = LocalElasticContactForce[0] * MaximumAdmisibleShearForce / ActualElasticShearForce;
                LocalElasticContactForce[1] = LocalElasticContactForce[1] * MaximumAdmisibleShearForce / ActualElasticShearForce;
                ViscoDampingLocalContactForce[0] = 0.0;
                ViscoDampingLocalContactForce[1] = 0.0;
            }
            
            else {
                
                const double ActualViscousShearForce = MaximumAdmisibleShearForce - ActualElasticShearForce;
                double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0]\
                                                        + ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);
        
                ViscoDampingLocalContactForce[0] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                ViscoDampingLocalContactForce[1] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
            }
                    
            sliding = true;    
        }
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall,
                                                                double indentation) {                                        
        
        const double my_mass               = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();
            
        double coefficient_of_restitution;
        
        if (my_ln_of_restit_coeff > 0.0) coefficient_of_restitution = 0.0;
        
        else coefficient_of_restitution = exp(my_ln_of_restit_coeff);
        
        //const double RCTS = 2.0 * sqrt(my_mass / mKn);
        
        //std::cout << "1% of Rayleigh's critical time step is = " << 0.01 * RCTS << std::endl;
        
        // NORMAL FORCE
        
        double normal_damping_coefficient;
        
        if (!my_ln_of_restit_coeff) normal_damping_coefficient = 0.0;
        
        else {
            //normal_damping_coefficient = (0.833333333333 * coefficient_of_restitution * coefficient_of_restitution - 2.25 * coefficient_of_restitution + 1.41666666666666) * sqrt(my_mass * mKn) * sqrt(sqrt(indentation));
            double gamma = 0.0;
            CalculateGamma(coefficient_of_restitution, gamma);        
            normal_damping_coefficient = 2.0 * gamma * sqrt(my_mass * mKn) * sqrt(sqrt(indentation));
        }
        
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient * LocalRelVel[2];
        
        // TANGENTIAL FORCE
        
        const double tangential_damping_coefficient = normal_damping_coefficient;
               
        //if (!sliding) { //tangential forces are not calculated yet, so sliding can never be true
        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        //}
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForce(const double indentation, SphericParticle* const element1, SphericParticle* const element2){
        return mKn* pow(indentation, 1.5);
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        return mKn* pow(indentation, 1.5);
    }

    double DEM_D_Tsuji_viscous_Coulomb::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2){        
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce(element1, element2);
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM(element, wall);
    }            
    
} /* namespace Kratos.*/
