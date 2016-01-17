//Authors: M.A. Celigueta and G. Casas (CIMNE)
//Date: January 2016

#include "DEM_D_Bentonite_Colloid_CL.h"

double ToThePower(double base, unsigned int n)
{
    if (n == 0) return 1.0;
    double x = base;

    while (n > 1){
        x *= base;
        n -= 1;
    }

    return  x;
}

double GetDebyeLength(double cation_concentration)
{
    return 1.04e8 / sqrt(cation_concentration);
}

namespace Kratos {
//G Hard-coded values for the moment, some should probably be nodal
    DEM_D_Bentonite_Colloid::DEM_D_Bentonite_Colloid(){
        mA_H = 10e-19;
        double d_p = 2.0e-7; // particle diameter; it whould be equal for both particles or the third law of Newton will be violated
        mA_p = 0.25 * KRATOS_M_PI * d_p * d_p;
        mThickness = 1.0e-9;
        mDDLCoefficient = 1.56e6;
        mDebyeLengthInv = 1.04e8;
        mEquivRadius = d_p / KRATOS_M_PI; // this is the "coin" equivalent radius
    }
//Z
    double mA_p;
    double mThickness;
    double mDDLCoefficient;

    void DEM_D_Bentonite_Colloid::Initialize(const ProcessInfo& rCurrentProcessInfo) {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Bentonite_Colloid::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Bentonite_Colloid(*this));
        return p_clone;
    }

    void DEM_D_Bentonite_Colloid::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEM_D_Bentonite_Colloid to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    std::string DEM_D_Bentonite_Colloid::GetTypeOfLaw() {
        std::string type_of_law = "Linear";
        return type_of_law;
    }

    ///////////////////////// 
    // DEM-DEM INTERACTION //
    /////////////////////////
    
    void DEM_D_Bentonite_Colloid::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
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
    
    void DEM_D_Bentonite_Colloid::CalculateForces(ProcessInfo& rCurrentProcessInfo,
                                                       const double OldLocalContactForce[3],
                                                             double LocalElasticContactForce[3],
                                                             double LocalDeltDisp[3],
                                                             double LocalRelVel[3],            
                                                             double indentation,
                                                             double previous_indentation,
                                                             double ViscoDampingLocalContactForce[3],
                                                             double& cohesive_force,
                                                             SphericParticle* element1,
                                                             SphericParticle* element2,
                                                             bool& sliding) {

        InitializeContact(element1, element2, indentation);
//G
        //LocalElasticContactForce[2]  = CalculateNormalForce(indentation);
        double distance = element1->GetRadius() + element2->GetRadius() - indentation;
        double cation_concentration = element1->GetGeometry()[0].FastGetSolutionStepValue(CATION_CONCENTRATION);//Z
        LocalElasticContactForce[2]  = CalculateNormalForce(distance, cation_concentration);
//Z
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2, indentation);                                                      
        
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        /*if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
	    }*/
        
        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element1, element2, indentation, previous_indentation);
        
    }
    
    void DEM_D_Bentonite_Colloid::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2) {
        
//        const double my_mass    = element1->GetMass();
//        const double other_mass = element2->GetMass();
                
//        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);
        
//        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
//        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
//        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
//        const double equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mKn);
//        const double equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mKt);
              
//        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
//        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
//        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];

        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;
    }
    
    ///////////////////////// 
    // DEM-FEM INTERACTION //
    /////////////////////////
    
    void DEM_D_Bentonite_Colloid::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta) {
                
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

        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius = sqrt(effective_radius);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }    
    
    void DEM_D_Bentonite_Colloid::CalculateForcesWithFEM(ProcessInfo& rCurrentProcessInfo,
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
//G
        double distance = element->GetRadius() - indentation;
        double cation_concentration = element->GetGeometry()[0].FastGetSolutionStepValue(CATION_CONCENTRATION);

        LocalElasticContactForce[2]  = CalculateNormalForce(distance, cation_concentration);
//Z
        cohesive_force              = CalculateCohesiveNormalForceWithFEM(element, wall, indentation);
        
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);
        
        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
                
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }
        
        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                        sliding, element, wall, indentation, previous_indentation);
    }
        
    template<class NeighbourClassType>
    void DEM_D_Bentonite_Colloid::CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                                      const double OldLocalContactForce[3],
                                                                      double LocalElasticContactForce[3],
                                                                      double ViscoDampingLocalContactForce[3],
                                                                      const double LocalDeltDisp[3],            
                                                                      bool& sliding,
                                                                      SphericParticle* const element,
                                                                      NeighbourClassType* const neighbour,
                                                                      double indentation,
                                                                      double previous_indentation) {
                                                                        
        LocalElasticContactForce[0] = 0.0;
        LocalElasticContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
    }
    
    void DEM_D_Bentonite_Colloid::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                SphericParticle* const element,
                                                                DEMWall* const wall) {                                        
        
        const double my_mass    = element->GetMass();              
        const double gamma = element->GetProperties()[DAMPING_GAMMA];
        const double normal_damping_coefficient     = 2.0 * gamma * sqrt(my_mass * mKn);
        const double tangential_damping_coefficient = 2.0 * gamma * sqrt(my_mass * mKt);
              
        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient     * LocalRelVel[2];  

    }
    
    double DEM_D_Bentonite_Colloid::CalculateNormalForce(const double distance, const double cation_concentration){
        
        double F_vdW = CalculateVanDerWaalsForce(distance);
        double F_DDL = CalculateDiffuseDoubleLayerForce(distance, cation_concentration);

        return F_vdW + F_DDL;
    }

    double DEM_D_Bentonite_Colloid::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation){
        return 0.0;
    }
    
    
    double DEM_D_Bentonite_Colloid::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation){
        return 0.0;
    }

    double DEM_D_Bentonite_Colloid::CalculateVanDerWaalsForce(const double distance)
    {
        return - mA_p * mA_H / (6.0 * KRATOS_M_PI) * (1.0 / ToThePower(distance, 3) - 2.0 / ToThePower(distance + mThickness, 3) + 1.0 / ToThePower(distance + 2 * mThickness, 3));
    }

    double DEM_D_Bentonite_Colloid::CalculateDiffuseDoubleLayerForce(const double distance, const double cation_concentration)
    {

        return mA_p * mDDLCoefficient * exp(- GetDebyeLength(cation_concentration) * distance);
    }

} // namespace Kratos
