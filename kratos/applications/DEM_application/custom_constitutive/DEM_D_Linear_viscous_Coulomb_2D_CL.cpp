// Authors: M.A. Celigueta and S. Latorre (CIMNE)
// Date: October 2015
// System includes
//#include <string>
//#include <iostream>

// External includes


// Project includes
//#include "DEM_application.h"
#include "DEM_D_Linear_viscous_Coulomb_2D_CL.h"
//#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Linear_viscous_Coulomb2D::Initialize(const ProcessInfo& rCurrentProcessInfo) {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_viscous_Coulomb2D::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_viscous_Coulomb2D(*this));
        return p_clone;
    }

    void DEM_D_Linear_viscous_Coulomb2D::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << " Assigning DEM_D_Linear_viscous_Coulomb2D to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    ///////////////////////// 
    // DEM-DEM INTERACTION //
    /////////////////////////
    
    void DEM_D_Linear_viscous_Coulomb2D::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        //Get equivalent Radius
        const double my_radius     = element1->GetRadius();
        const double other_radius  = element2->GetRadius();
        const double radius_sum    = my_radius + other_radius;
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
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);                   
        
        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius = sqrt(equiv_radius);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;       
        
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateForces(ProcessInfo& rCurrentProcessInfo,
                                                       const double OldLocalContactForce[3],
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
        
        InitializeContact(element1, element2, indentation);
        
        LocalElasticContactForce[2]  = CalculateNormalForce(indentation);
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
                double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] + ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);
        
                ViscoDampingLocalContactForce[0] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                ViscoDampingLocalContactForce[1] *= ActualViscousShearForce / ViscoDampingLocalContactForceModule;
            }
                    
            sliding = true;   
        }               
    }
    void DEM_D_Linear_viscous_Coulomb2D::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2) {
        
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
    
    double DEM_D_Linear_viscous_Coulomb2D::CalculateNormalForce(const double indentation) {
        //return 0.666666666666666666667 * mKn * indentation;
        return mKn* indentation;
    }

    double DEM_D_Linear_viscous_Coulomb2D::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation){        
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce(element1, element2, indentation);
    }
    
    
    ///////////////////////// 
    // DEM-FEM INTERACTION //
    /////////////////////////
    
    void DEM_D_Linear_viscous_Coulomb2D::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta) {
                
        const double my_radius        = element->GetRadius(); // Get equivalent Radius
        const double effective_radius = my_radius - ini_delta;
        const double my_young         = element->GetYoung(); // Get equivalent Young's Modulus
        const double walls_young      = wall->GetYoung();
        const double my_poisson       = element->GetPoisson();
        const double walls_poisson    = wall->GetPoisson();
                
        const double equiv_young      = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));
                
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus); 
        /*
        const double effective_young = my_young / (1.0 - my_poisson * my_poisson); // Equivalent Young Modulus for RIGID WALLS! 
        const double effective_young = 0.5 * my_young / (1.0 - my_poisson * my_poisson); // Equivalent Young Modulus if the wall has the same E of the sphere
        // Infinite Young was imposed to the wall in the formula that includes both. 
        // For flexible walls, the computation is different, Thornton 2012 proposes same Young for both)
        
        const double effective_shear = 0.5 * my_young / ((1.0 + my_poisson) * (2.0 - my_poisson)); // Equivalent Shear Modulus for RIGID WALLS!
        const double effective_shear = 0.25 * my_young / ((1.0 + my_poisson) * (2.0 - my_poisson)); // Equivalent Shear Modulus if the wall has the same E of the sphere
        // Infinite Shear was imposed to the wall in the formula that includes both. 
        // For flexible walls, the computation is different, Thornton 2012 proposes same Shear for both)
        */        
        //Normal and Tangent elastic constants
        //const double sqrt_equiv_radius_and_indentation = sqrt(effective_radius * indentation);
        //mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        const double sqrt_equiv_radius = sqrt(effective_radius);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }    
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateForcesWithFEM(ProcessInfo& rCurrentProcessInfo,
                                                              const double OldLocalContactForce[3],
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
        InitializeContactWithFEM(element, wall, indentation);
        
        LocalElasticContactForce[2] = CalculateNormalForce(indentation);
        cohesive_force              = CalculateCohesiveNormalForceWithFEM(element, wall, indentation);                                                      
        
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, sliding, element, wall, indentation);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForceWithFEM(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element, wall, indentation, previous_indentation);
        
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

        const double my_tg_of_friction_angle    = element->GetTgOfFrictionAngle();
        const double wall_tg_of_friction_angle  = wall->GetTgOfFrictionAngle();
        const double equiv_tg_of_fri_ang        = 0.5 * (my_tg_of_friction_angle + wall_tg_of_friction_angle);    
        
        const double MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_fri_ang;
        
        const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];
        
        const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);
                
        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {
            
            const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            
            const double dot_product = LocalElasticContactForce[0] * ViscoDampingLocalContactForce[0] + LocalElasticContactForce[1] * ViscoDampingLocalContactForce[1];
            const double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +\
                                                                    ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);
            if (dot_product >= 0.0) {
                if (ActualElasticShearForce > MaximumAdmisibleShearForce) {
                    const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                    LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                    LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                ViscoDampingLocalContactForce[0] = 0.0;
                ViscoDampingLocalContactForce[1] = 0.0;
            }
            
            else {
                
                const double ActualViscousShearForce = MaximumAdmisibleShearForce - ActualElasticShearForce;
                    const double fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                    ViscoDampingLocalContactForce[0]    *= fraction;
                    ViscoDampingLocalContactForce[1]    *= fraction;
            }
                    
            }
            else {
                if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                    const double fraction = (MaximumAdmisibleShearForce + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                    ViscoDampingLocalContactForce[0] *= fraction;
                    ViscoDampingLocalContactForce[1] *= fraction;
                }
                else {
                    const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                    LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                    LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                    ViscoDampingLocalContactForce[0] = 0.0;
                    ViscoDampingLocalContactForce[1] = 0.0;
                }
            }
            sliding = true;    
        }
    }
    
    void DEM_D_Linear_viscous_Coulomb2D::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall,
                                                                double indentation){                                        
        
        const double my_mass    = element->GetMass();              
        const double gamma = element->GetProperties()[DAMPING_GAMMA];
        const double normal_damping_coefficient     = 2.0 * gamma * sqrt(my_mass * mKn);
        const double tangential_damping_coefficient = 2.0 * gamma * sqrt(my_mass * mKt);
              
        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient     * LocalRelVel[2];  

    }
    
    double DEM_D_Linear_viscous_Coulomb2D::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation){
        return DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM(element, wall, indentation);
    }

} /* namespace Kratos.*/
