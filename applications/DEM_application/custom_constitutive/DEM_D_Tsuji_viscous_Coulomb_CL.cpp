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
    
    
    void DEM_D_Tsuji_viscous_Coulomb::InitializeContact(SphericParticle* const element1,
                                                 SphericParticle* const element2) {
        //Get equivalent Radius
        const double my_radius     = element1->GetRadius();
        const double other_radius  = element2->GetRadius();
        const double radius_sum    = my_radius + other_radius;
        const double radius_sum_i  = 1 / radius_sum;
        const double equiv_radius  = my_radius * other_radius * radius_sum_i;
        
        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        
        //Get equivalent Poisson's Modulus
        double equiv_poisson;
        if ((my_poisson + other_poisson)!= 0.0) {
            equiv_poisson            = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson);
        } else {
            equiv_poisson            = 0.0;
        }                  
        
        //Normal and Tangent elastic constants
        mKn = 1.3333333333333333 * equiv_young * sqrt(equiv_radius);
        mKt = 2.0 * mKn * (1.0 - equiv_poisson * equiv_poisson) / ((2.0 - equiv_poisson) * (1.0 + equiv_poisson));
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
                                                                      const double ViscoDampingLocalContactForce[3],
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
        const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        
        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {

            LocalElasticContactForce[0] = LocalElasticContactForce[0] * MaximumAdmisibleShearForce / ActualElasticShearForce;
            LocalElasticContactForce[1] = LocalElasticContactForce[1] * MaximumAdmisibleShearForce / ActualElasticShearForce;

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
        
        if (my_ln_of_restit_coeff == 1.0) coefficient_of_restitution = 0.0;
        
        else coefficient_of_restitution = exp(my_ln_of_restit_coeff);
        
        const double RCTS = 2.0 * sqrt(my_mass / mKn);
        
        std::cout << "1% of Rayleigh's critical time step is = " << 0.01 * RCTS << std::endl;
        
        // NORMAL FORCE
        
        double normal_damping_coefficient;
        
        if (!my_ln_of_restit_coeff) normal_damping_coefficient = 0.0;
        
        else normal_damping_coefficient = (0.833333333333 * coefficient_of_restitution * coefficient_of_restitution - 2.25 * coefficient_of_restitution + 1.41666666666666) * sqrt(my_mass * mKn) * sqrt(sqrt(indentation));
        
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient * LocalRelVel[2];
        
        // TANGENTIAL FORCE
        
        const double tangential_damping_coefficient = normal_damping_coefficient;
               
        if (!sliding) {
            ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        }
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForce(const double indentation, SphericParticle* const element1, SphericParticle* const element2){
        return mKn* pow(indentation, 1.5);
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        return mKn* pow(indentation, 1.5);
    }

    double DEM_D_Tsuji_viscous_Coulomb::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2){        
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce( element1, element2);
    }
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM( element, wall);
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateTangentialForce(const double normal_force,
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element1,
                                                    SphericParticle* const element2) {                                
        DEMDiscontinuumConstitutiveLaw::CalculateStandardTangentialForce(normal_force, LocalElasticContactForce, LocalDeltDisp,  sliding, element1, element2);        
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2){                                        
        
        DEMDiscontinuumConstitutiveLaw::CalculateStandardViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, sliding, element1, element2);  
    }
    
} /* namespace Kratos.*/
