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
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateForces(double LocalElasticContactForce[3],
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
        
        CalculateNormalForceTsuji(LocalElasticContactForce, kn_el, indentation);
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

    void DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForceTsuji(double LocalElasticContactForce[3], double kn_el, double indentation) {
        LocalElasticContactForce[2] = kn_el * pow(indentation, 1.5); //  1 --- Tsujiian (non-linear compression, linear tension)
    }

    void DEM_D_Tsuji_viscous_Coulomb::CalculateTangentialForceLinear(double LocalElasticContactForce[3],
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

    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDamping(double LocalRelVel[3],
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
        mKn              = 1.3333333333333333 * equiv_young * sqrt(equiv_radius);
        mKt              = 2.0 * mKn * (1.0 - equiv_poisson * equiv_poisson) / ((2.0 - equiv_poisson) * (1.0 + equiv_poisson));
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta) {
        //Get equivalent Radius
        const double my_radius     = element->GetRadius();
        
        //Get equivalent Young's Modulus
        const double my_young        = element->GetYoung();
        const double my_poisson      = element->GetPoisson();
                
        //Normal and Tangent elastic constants
        mKn = (1.33333333 * my_young / (1.0 - my_poisson * my_poisson)) * sqrt(my_radius - ini_delta); //TODO: CHECK THIS 
        mKt = 8.0 * sqrt(my_radius) * my_young * 0.5 / ((1.0 + my_poisson) * (2.0 - my_poisson)); //LATER ON THIS VALUE WILL BE MULTIPLIED BY SQRT(INDENTATION)                       
        
    }
    
    
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForce(const double indentation, SphericParticle* const element1, SphericParticle* const element2){
        //KRATOS_WATCH(mKn* pow(indentation, 1.5))
        return mKn* pow(indentation, 1.5);
    }
    double DEM_D_Tsuji_viscous_Coulomb::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        //KRATOS_WATCH(mKn* pow(indentation, 1.5))
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
    void DEM_D_Tsuji_viscous_Coulomb::CalculateTangentialForceWithFEM(const double normal_force,
                                                                    double LocalElasticContactForce[3],
                                                                    const double LocalDeltDisp[3],            
                                                                    bool& sliding,
                                                                    SphericParticle* const element,
                                                                    DEMWall* const wall,
                                                                    double indentation) {
                                                                        
        const double WallBallFriction = wall->mTgOfFrictionAngle; //TODO: CHECK THIS (SHOULD IT BE AN AVERAGE?)      
        
        LocalElasticContactForce[0] += - mKt * sqrt(indentation) * LocalDeltDisp[0];
        LocalElasticContactForce[1] += - mKt * sqrt(indentation) * LocalDeltDisp[1];

        const double ShearForceMax = normal_force * WallBallFriction;
        const double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        if (ShearForceNow > ShearForceMax) {
            double inv_ShearForceNow = 1.0 / ShearForceNow;
            LocalElasticContactForce[0] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[1];
            sliding = true;
        }                                                                                                
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2){                                        
        
        DEMDiscontinuumConstitutiveLaw::CalculateStandardViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, sliding, element1, element2);  
    }
    
    void DEM_D_Tsuji_viscous_Coulomb::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall,
                                                                double indentation) {                                        
        
        const double my_mass               = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();
                
        double e = exp(my_ln_of_restit_coeff);
        
        // NORMAL FORCE
        
        double nu_normal = 1.0 * (0.833333333333 * e * e - 2.25 * e + 1.41666666666666) * sqrt(my_mass * mKn) * sqrt(sqrt(indentation));
        
        ViscoDampingLocalContactForce[2] = - nu_normal * LocalRelVel[2];
        
        // TANGENTIAL FORCE
        
        double nu_tangential = nu_normal;
               
        if (!sliding) {
            ViscoDampingLocalContactForce[0] = - nu_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = - nu_tangential * LocalRelVel[1];
        }
    }
    
} /* namespace Kratos.*/
