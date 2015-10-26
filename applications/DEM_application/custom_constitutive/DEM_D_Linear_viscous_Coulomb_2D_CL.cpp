// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Linear_viscous_Coulomb2D::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_viscous_Coulomb2D::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_viscous_Coulomb2D(*this));
        return p_clone;
    }

    void DEM_D_Linear_viscous_Coulomb2D::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << " Assigning DEM_D_Linear_viscous_Coulomb2D to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateContactArea(double radius, double other_radius, double &calculation_area) {
        
        KRATOS_TRY
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = KRATOS_M_PI * rmin*rmin;
        KRATOS_CATCH("")  
    }

    void DEM_D_Linear_viscous_Coulomb2D::CalculateElasticConstants(double &kn_el,
            double &kt_el,
            double initial_dist,
            double equiv_young,
            double equiv_poisson,
            double calculation_area) {
        
        KRATOS_TRY 
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")  
    }
    
    
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateForces(double LocalElasticContactForce[3],
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
        CalculateNormalForceLinear(LocalElasticContactForce, kn_el, indentation);
        //CalculateCohesiveForce(cohesion_force, cohesion, cohesion_area);
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

    void DEM_D_Linear_viscous_Coulomb2D::CalculateNormalForceLinear(double LocalElasticContactForce[3], const double kn_el, const double indentation) {
        
        KRATOS_TRY
        
        LocalElasticContactForce[2] = kn_el * indentation; //  0 ---linear compression & tension
        KRATOS_CATCH("")        
    }

    void DEM_D_Linear_viscous_Coulomb2D::CalculateTangentialForceLinear(double LocalElasticContactForce[3],
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

    void DEM_D_Linear_viscous_Coulomb2D::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding,
            int mDampType) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in oposite direction the visco damping can't overpass the force...
        
        KRATOS_TRY
                
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
    
    void DEM_D_Linear_viscous_Coulomb2D::InitializeContact(SphericParticle* const element1,
                                                 SphericParticle* const element2) {
        //Get equivalent Radius
        const double my_radius     = element1->GetRadius();
        const double other_radius  = element2->GetRadius();
        const double radius_sum    = my_radius + other_radius;
        const double radius_sum_i  = 1 / radius_sum;
        const double equiv_radius  = 2 * my_radius * other_radius * radius_sum_i;
        
        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young = 2.0 * my_young * other_young / (my_young + other_young);
        
        //Get equivalent Poisson's Modulus
        double equiv_poisson;
        if ((my_poisson + other_poisson)!= 0.0) {
            equiv_poisson            = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson);
        } else {
            equiv_poisson            = 0.0;
        }                  
        
        double equiv_area = 0.25 * KRATOS_M_PI * equiv_radius * equiv_radius; // 0.25 is because we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.        
        
        //Normal and Tangent elastic constants
        mKn              = equiv_young * equiv_area * radius_sum_i; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
        
        //mKt              = mKn / (2.0 + equiv_poisson + equiv_poisson);                               
        mKt = mKn * 2 * (1-equiv_poisson) / (2-equiv_poisson); //Cleary, Thornton, based on Hertzian Theory
        
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateForces(const double OldLocalContactForce[3],
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
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2, indentation);                                                      
        
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForce(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element1, element2, indentation, previous_indentation);
        
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta) {
        //Get equivalent Radius
        const double my_radius     = element->GetRadius();
        
        //Get equivalent Young's Modulus
        const double my_young        = element->GetYoung();
        const double my_poisson      = element->GetPoisson();
        const double area            = KRATOS_M_PI * my_radius * my_radius;                       
                
        //Normal and Tangent elastic constants
        mKn = my_young * area / (2.0 * my_radius - ini_delta);        
        mKt = mKn / (2.0 * (1.0 + my_poisson));                                
    }    
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateForcesWithFEM(const double OldLocalContactForce[3],
                                                             double LocalElasticContactForce[3],
                                                             double LocalDeltDisp[3],
                                                             double LocalRelVel[3],            
                                                             double indentation,
                                                             double previous_indentation,
                                                             double ViscoDampingLocalContactForce[3],
                                                             double& cohesive_force,
                                                             SphericParticle* const element,
                                                             DEMWall* const wall,
                                                             bool& sliding) {
        InitializeContactWithFEM(element, wall);
        
        LocalElasticContactForce[2]  = CalculateNormalForceWithFEM(indentation, element, wall);
        cohesive_force               = CalculateCohesiveNormalForceWithFEM(element, wall, indentation);                                                      
        
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, sliding, element, wall, indentation);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForceWithFEM(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element, wall, indentation, previous_indentation);
        
        //if (sliding) ViscoDampingLocalContactForce[0] = ViscoDampingLocalContactForce[1] = 0.0;
    }
    
    double DEM_D_Linear_viscous_Coulomb2D::CalculateNormalForce(const double indentation, SphericParticle* const element1, SphericParticle* const element2){
        return mKn* indentation;
    }
    double DEM_D_Linear_viscous_Coulomb2D::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        return mKn* indentation;
    }

    double DEM_D_Linear_viscous_Coulomb2D::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce(element1, element2, indentation);                
    }
    double DEM_D_Linear_viscous_Coulomb2D::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM( element, wall, indentation);
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateTangentialForce(const double normal_contact_force,
                                                                const double OldLocalContactForce[3],
                                                                double LocalElasticContactForce[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                const double LocalDeltDisp[3],            
                                                                bool& sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2,
                                                                double indentation,
                                                                double previous_indentation) { 
        
        
        LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * LocalDeltDisp[1];

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
    void DEM_D_Linear_viscous_Coulomb2D::CalculateTangentialForceWithFEM(const double normal_contact_force,
                                                                      const double OldLocalContactForce[3],
                                                                      double LocalElasticContactForce[3],
                                                                      double ViscoDampingLocalContactForce[3],
                                                                      const double LocalDeltDisp[3],            
                                                                      bool& sliding,
                                                                      SphericParticle* const element,
                                                                      DEMWall* const wall,
                                                                      double indentation,
                                                                      double previous_indentation) {
                                                                        
        LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * LocalDeltDisp[1];
        
        const double WallBallFriction = wall->GetTgOfFrictionAngle();
        
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
    void DEM_D_Linear_viscous_Coulomb2D::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2){                                        
        
        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
                
        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);        
        
        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
        const double equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mKn);
        const double equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mKt);                
              
        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];  
               
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall,
                                                                double indentation){                                        
        
        const double my_mass    = element->GetMass();              
        const double gamma = element->GetProperties()[DAMPING_GAMMA];
        const double equiv_visco_damp_coeff_normal     = 2.0 * gamma * sqrt(my_mass * mKn);
        const double equiv_visco_damp_coeff_tangential = 2.0 * gamma * sqrt(my_mass * mKt);                
              
        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];  

    }
    

} /* namespace Kratos.*/
