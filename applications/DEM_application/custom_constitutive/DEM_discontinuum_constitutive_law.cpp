// Last modified by: S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_application.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw() {
        //std::cout << " DEMDiscontinuumConstitutiveLaw constructor..." << std::endl;

    } // Class DEMDiscontinuumConstitutiveLaw

    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw) {
        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    }

    //    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw){
    //        //std::cout << " DEMDiscontinuumConstitutiveLaw copy constructor..." << std::endl;
    //    }

    void DEMDiscontinuumConstitutiveLaw::Initialize(const ProcessInfo& r_process_info) {
    }

    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        //std::cout << "Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }
    
    std::string DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw() {
        //std::cout << "Law destructor..." ; 
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateContactArea(double radius, double other_radius, double &calculation_area) {
        
        KRATOS_TRY
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = KRATOS_M_PI * rmin*rmin;
        KRATOS_CATCH("")  
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateElasticConstants(double &kn_el,
            double &kt_el,
            double initial_dist,
            double equiv_young,
            double equiv_poisson,
            double calculation_area) {
        
        KRATOS_TRY 
        double equiv_shear = equiv_young / (2.0 * (1.0 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")  
    }    
    
    void DEMDiscontinuumConstitutiveLaw::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta) {        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::InitializeContact) should not be called.","")
    }
    
    void DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM) should not be called.","")
    }
    
    void DEMDiscontinuumConstitutiveLaw::GetContactStiffness(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta, double& kn,double& kt){
      
      InitializeContact(element1, element2, ini_delta);
      kn = mKn;
      kt = mKt;
      
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateForces(ProcessInfo& r_process_info,
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
            
        LocalElasticContactForce[2]  = CalculateNormalForce(indentation, element1, element2);   
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2, indentation);

        CalculateTangentialForce(LocalElasticContactForce[2], LocalElasticContactForce, LocalDeltDisp, sliding, element1, element2);        
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                SphericParticle* element1,
                                                                SphericParticle* element2) {  
        
        
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateForcesWithFEM(ProcessInfo& r_process_info,
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
                
        LocalElasticContactForce[2]  = CalculateNormalForceWithFEM(indentation, element, wall);
        cohesive_force               = CalculateCohesiveNormalForceWithFEM(element, wall, indentation);                                                      

        CalculateTangentialForceWithFEM(OldLocalContactForce, LocalElasticContactForce, LocalDeltDisp, sliding, element, wall, indentation, previous_indentation);               
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(const double indentation, SphericParticle* element1, SphericParticle* element2) {        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) should not be called.","")
        return 0.0;
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(const double indentation) {        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) should not be called.","")
        return 0.0;
    }
        
    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForceWithFEM) should not be called.","")
        return 0.0;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateTangentialForce(const double normal_force,
                                                                    double LocalElasticContactForce[3],
                                                                    const double LocalDeltDisp[3],            
                                                                    bool& sliding,
                                                                    SphericParticle* const element1,
                                                                    SphericParticle* const element2) {        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateTangentialForce) should not be called.","")
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateTangentialForceWithFEM(const double OldLocalContactForce[3],
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation,
                                                    double previous_indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateTangentialForceWithFEM) should not be called.","")
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {        
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce) should not be called.","")
        return 0.0;        
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation){
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM) should not be called.","")
        return 0.0;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                    double ViscoDampingLocalContactForce[3],
                                                                    SphericParticle* element1,
                                                                    SphericParticle* element2) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForce) should not be called.","")
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForceWithFEM(double LocalRelVel[3], double ViscoDampingLocalContactForce[3], SphericParticle* const element, DEMWall* const wall){
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForceWithFEM) should not be called.","")
    }
    
    
    ////////////////////////////////////////////////
    /////// LIBRARY OF CALCULATION FUNCTIONS ///////
    ////////////////////////////////////////////////                
              
    void DEMDiscontinuumConstitutiveLaw::CalculateForces(ProcessInfo& r_process_info,
                                                         double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double kn_el,
            double kt_el,
            double indentation,
            double& failure_criterion_state,
            bool& sliding,
            SphericParticle* element1,
            SphericParticle* element2,
            int &mNeighbourFailureId_count) {

        KRATOS_TRY
        CalculateNormalForceLinear(LocalElasticContactForce, kn_el, indentation);
        CalculateTangentialForceLinear(LocalElasticContactForce,
                LocalDeltDisp,
                kt_el,
                indentation,
                failure_criterion_state,
                sliding,
                element1,
                element2,
                mNeighbourFailureId_count);
        KRATOS_CATCH("")  
        
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateNormalForceLinear(double LocalElasticContactForce[3], const double kn_el, const double indentation) {
        KRATOS_TRY
        LocalElasticContactForce[2] = kn_el * indentation; //  0 ---linear compression & tension
        KRATOS_CATCH("")  
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateNormalForceHertz(double LocalElasticContactForce[3], double kn_el, double indentation) {
        KRATOS_TRY
        LocalElasticContactForce[2] = kn_el * pow(indentation, 1.5); //  1 --- Hertzian (non-linear compression, linear tension)
        KRATOS_CATCH("")  
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateTangentialForceLinear(double LocalElasticContactForce[3],
            double LocalDeltDisp[3],
            double kt_el,
            double indentation,
            double& failure_criterion_state,
            bool& sliding,
            SphericParticle* element1,
            SphericParticle* element2,
            int &mNeighbourFailureId_count) {

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

    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in opposite direction the visco damping can't overpass the force...

        KRATOS_TRY  

                ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];


            if (sliding == false) { //only applied when no sliding to help to the regularized friction law or the spring convergence
                ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
                ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
            }

    KRATOS_CATCH("")      
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            SphericParticle* element1,
            SphericParticle* element2,
            double kn_el,
            double kt_el) {
        
        KRATOS_TRY 
        double aux_norm_to_tang = 0.0;
        const double my_mass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
//        const double mDempack_local_damping = element1->GetProperties()[DEMPACK_LOCAL_DAMPING];
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];

        equiv_visco_damp_coeff_normal = (1-mCoefficientOfRestitution) * 2.0 * sqrt(kn_el / (my_mass + other_real_mass)) * (sqrt(my_mass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
        KRATOS_CATCH("")  
    }   
    
} // KRATOS
