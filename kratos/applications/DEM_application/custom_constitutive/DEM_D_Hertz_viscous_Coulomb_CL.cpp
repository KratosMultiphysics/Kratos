// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Hertz_viscous_Coulomb::Initialize(const ProcessInfo& rCurrentProcessInfo){
        
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Hertz_viscous_Coulomb::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Hertz_viscous_Coulomb(*this));
        return p_clone;
    }

    void DEM_D_Hertz_viscous_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        //std::cout << " Assigning DEM_D_Hertz_viscous_Coulomb to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateForces(double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double kn_el,
            double kt_el,
            double indentation,
            double& failure_criterion_state,
            bool& sliding,
            SphericParticle* element1,
            SphericParticle* element2,
            int &mNeighbourFailureId_count,
            double mapping_new_cont) {
        
        KRATOS_TRY
        
        CalculateNormalForceHertz(LocalElasticContactForce, kn_el, indentation);
        CalculateTangentialForceLinear(LocalElasticContactForce,
                LocalDeltDisp,
                kt_el,
                indentation,
                failure_criterion_state,
                sliding,
                element1,
                element2,
                mNeighbourFailureId_count,
                mapping_new_cont);    
        
        KRATOS_CATCH("")
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateNormalForceHertz(double LocalElasticContactForce[3], double kn_el, double indentation) {
        LocalElasticContactForce[2] = kn_el * pow(indentation, 1.5); 
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateTangentialForceLinear(double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double kt_el,
            double indentation,
            double& failure_criterion_state,
            bool& sliding,
            SphericParticle* element1,
            SphericParticle* element2,
            int &mNeighbourFailureId_count,
            double mapping_new_cont) {
        
        
        KRATOS_TRY
        const double other_tg_of_fri_angle = element2->GetTgOfFrictionAngle();
        const double myTgOfFrictionAngle = element1->GetTgOfFrictionAngle();

        if (mNeighbourFailureId_count != 0) //*   //degut als canvis de DEMPACK hi ha hagut una modificació, ara despres de trencar es fa akest maping de maxima tangencial que és correcte!
        {
            LocalElasticContactForce[0] += -1.0 * kt_el * LocalDeltDisp[0]; // 0: first tangential
            LocalElasticContactForce[1] += -1.0 * kt_el * LocalDeltDisp[1]; // 1: second tangential  

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            double Frictional_ShearForceMax = (0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle)) * LocalElasticContactForce[2];
            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            failure_criterion_state = 1.0;
            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }
        KRATOS_CATCH("")
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding,
            int mDampType) {
        
        KRATOS_TRY

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in oposite direction the visco damping can't overpass the force...

        if (mDampType > 0) {

            if (mDampType == 11 || mDampType == 10) {
                ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];
            }

            if (sliding == false && (mDampType == 1 || mDampType == 11)) { //only applied when no sliding to help to the regularized friction law or the spring convergence
                ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
                ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
            }
        }
             KRATOS_CATCH("")   
    }
        
    
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////////NEW IMPLEMENTATION/////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    
    
    void DEM_D_Hertz_viscous_Coulomb::InitializeContact(SphericParticle* const element1,
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
    
    void DEM_D_Hertz_viscous_Coulomb::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta) {
        //Get equivalent Radius
        const double my_radius     = element->GetRadius();
        
        //Get equivalent Young's Modulus
        const double my_young           = element->GetYoung();
        const double my_poisson         = element->GetPoisson();
        const double my_shear_mod       = my_young / (2.0 * (1.0 + my_poisson));
        const double wall_young         = my_young;
        const double wall_poisson       = my_poisson;
        const double wall_shear_mod       = wall_young / (2.0 * (1.0 + wall_poisson));
        
        const double equiv_young              = my_young * wall_young / (wall_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - wall_poisson * wall_poisson));
        
        const double equiv_shear_mod = my_shear_mod * wall_shear_mod / (wall_shear_mod * (2.0 - my_poisson) + my_shear_mod * (2.0 - wall_poisson));
        
        //Normal and Tangent elastic constants
        mKn = 1.3333333333333 * equiv_young * sqrt(my_radius - ini_delta); //TODO: CHECK THIS 
        mKt = 8.0 * equiv_shear_mod * sqrt(my_radius * indentation);                       
        
    }
    
    void DEM_D_Hertz_viscous_Coulomb::CalculateForcesWithFEM(const double OldLocalContactForce[3],
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
        
        InitializeContactWithFEM(element, wall, indentation);
        
        double delta_force_normal = 0.0;
        
        LocalElasticContactForce[2]  = CalculateNormalForceWithFEM(indentation, element, wall);
        cohesive_force               = CalculateCohesiveNormalForceWithFEM(element, wall);                                                      
        
        const double my_mass               = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();
        const double my_radius             = element->GetRadius();
        const double my_young              = element->GetYoung();
        const double my_poisson            = element->GetPoisson();
        
        const double gamma = - my_ln_of_restit_coeff / sqrt(KRATOS_M_PI * KRATOS_M_PI + my_ln_of_restit_coeff * my_ln_of_restit_coeff);
        
        const double wall_young         = my_young;
        const double wall_poisson       = my_poisson;
        const double equiv_young        = my_young * wall_young / (wall_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - wall_poisson * wall_poisson));
        
        const double Kn = 2.0 * equiv_young * sqrt(my_radius * indentation);
        
        //double RCTS = 2.0 * sqrt(my_mass / mKn);
        
        //std::cout << "Rayleigh's critical time step is = " << RCTS << std::endl << "1% of Rayleigh's critical time step is = " << 0.01 * RCTS << std::endl;
        
        ViscoDampingLocalContactForce[2] = -2.0 * sqrt(my_mass * Kn) * gamma * LocalRelVel[2]; ///////////////////////////// vel or 2*vel??????
        
        delta_force_normal = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2] - OldLocalContactForce[2];
        
        CalculateTangentialForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, OldLocalContactForce, LocalElasticContactForce,
                                        LocalDeltDisp, sliding, element, wall, indentation, previous_indentation, delta_force_normal);
        
    }
    
    double DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(const double indentation, SphericParticle* const element1, SphericParticle* const element2){
        return mKn* pow(indentation, 1.5);
    }
    
    double DEM_D_Hertz_viscous_Coulomb::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        return mKn* pow(indentation, 1.5);
    }

    double DEM_D_Hertz_viscous_Coulomb::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2){        
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce( element1, element2);
    }
    
    double DEM_D_Hertz_viscous_Coulomb::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM( element, wall);
    }
    
    void DEM_D_Hertz_viscous_Coulomb::CalculateTangentialForce(const double normal_force,
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element1,
                                                    SphericParticle* const element2) {                                
        DEMDiscontinuumConstitutiveLaw::CalculateStandardTangentialForce(normal_force, LocalElasticContactForce, LocalDeltDisp,  sliding, element1, element2);        
    }
    
    void DEM_D_Hertz_viscous_Coulomb::CalculateTangentialForceWithFEM(double LocalRelVel[3],
                                                    double ViscoDampingLocalContactForce[3],
                                                    const double OldLocalContactForce[3],
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation,
                                                    double previous_indentation,
                                                    double delta_force_normal) {
        
        const double WallBallFriction = wall->mTgOfFrictionAngle;         //TODO: CHECK THIS (SHOULD IT BE AN AVERAGE?)      
                
        if (delta_force_normal < 0.0) {
            const double kt_quotient = sqrt(indentation / previous_indentation);
            LocalElasticContactForce[0] = LocalElasticContactForce[0] * kt_quotient - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = LocalElasticContactForce[1] * kt_quotient - mKt * LocalDeltDisp[1];
        }
        
        else {
            LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * LocalDeltDisp[1];
        }

        const double ShearForceMax = (LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2])* WallBallFriction;
        
        const double my_mass               = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();
        
        const double gamma = -1.0 * my_ln_of_restit_coeff / (sqrt(KRATOS_M_PI * KRATOS_M_PI + my_ln_of_restit_coeff * my_ln_of_restit_coeff));
        
        ViscoDampingLocalContactForce[0] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[1];
        
        double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];
        
        const double ShearForceNow = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

        if (ShearForceNow > ShearForceMax) {
            double normalizer = ShearForceMax / ShearForceNow;
            LocalElasticContactForce[0] = normalizer * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = normalizer * LocalElasticContactForce[1];
            sliding = true;
        }
        
        //if (!sliding ) { 
        //    ViscoDampingLocalContactForce[0] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[0];
        //    ViscoDampingLocalContactForce[1] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[1];
        //}
    }
    
    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2){                                        
        
        DEMDiscontinuumConstitutiveLaw::CalculateStandardViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, sliding, element1, element2);  
    }
    
    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall,
                                                                double indentation) {                                        
        
        const double my_mass               = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();
        
        const double gamma = -1.0 * my_ln_of_restit_coeff / (KRATOS_M_PI * KRATOS_M_PI + my_ln_of_restit_coeff * my_ln_of_restit_coeff);
        
        if (!sliding ) { 
            ViscoDampingLocalContactForce[0] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = - 2.0 * gamma * sqrt(my_mass * mKt) * LocalRelVel[1];
        }
    }
    
} /* namespace Kratos.*/
