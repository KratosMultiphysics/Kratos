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

    void DEMDiscontinuumConstitutiveLaw::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    }

    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        //std::cout << " Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw() {
        //std::cout << "Law destructor..." ; 
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateContactArea(double mRadius, double other_radius, double &calculation_area) {
        
        KRATOS_TRY
        double rmin = mRadius;
        if (other_radius < mRadius) rmin = other_radius;
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
    
    void DEMDiscontinuumConstitutiveLaw::InitializeContact(SphericParticle* element1, SphericParticle* element2) {        
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::InitializeContact) should not be called."<<std::endl<<std::flush;
    }
    
    void DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta){
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM) should not be called."<<std::endl<<std::flush;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateForces(double LocalElasticContactForce[3],
                                                        double LocalDeltDisp[3],
                                                        double LocalRelVel[3],            
                                                        double indentation,
                                                        double ViscoDampingLocalContactForce[3],
                                                        double& cohesive_force,
                                                        SphericParticle* element1,
                                                        SphericParticle* element2) {  
        
        bool sliding = false;
        InitializeContact(element1, element2);   
            
        LocalElasticContactForce[2]  = CalculateNormalForce(indentation, element1, element2);   
        cohesive_force               = CalculateCohesiveNormalForce(element1, element2);

        CalculateTangentialForce(LocalElasticContactForce[2], LocalElasticContactForce, LocalDeltDisp, sliding, element1, element2);        
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, sliding, element1, element2);   
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateForcesWithFEM(const double OldLocalContactForce[3],
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

        CalculateTangentialForceWithFEM(OldLocalContactForce, LocalElasticContactForce, LocalDeltDisp, sliding, element, wall, indentation, previous_indentation);               
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, sliding, element, wall, indentation);
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(double indentation,  SphericParticle* element1, SphericParticle* element2) {        
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) should not be called."<<std::endl<<std::flush;
        return 0.0;
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall){
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForceWithFEM) should not be called."<<std::endl<<std::flush;
        return 0.0;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateTangentialForce(const double normal_force,
                                                                    double LocalElasticContactForce[3],
                                                                    const double LocalDeltDisp[3],            
                                                                    bool& sliding,
                                                                    SphericParticle* const element1,
                                                                    SphericParticle* const element2) {        
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateTangentialForce) should not be called."<<std::endl<<std::flush;        
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateTangentialForceWithFEM(const double OldLocalContactForce[3],
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation,
                                                    double previous_indentation) {
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateTangentialForceWithFEM) should not be called."<<std::endl<<std::flush;        
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2) {        
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce) should not be called."<<std::endl<<std::flush;
        return 0.0;        
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall){
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM) should not be called."<<std::endl<<std::flush;
        return 0.0;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForce(double LocalRelVel[3],
                                                                    double ViscoDampingLocalContactForce[3],
                                                                    bool sliding,
                                                                    SphericParticle* element1,
                                                                    SphericParticle* element2) {
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForce) should not be called."<<std::endl<<std::flush;        
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                    double ViscoDampingLocalContactForce[3],
                                                    bool sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation){
        std::cout<<"This function (DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingForceWithFEM) should not be called."<<std::endl<<std::flush;        
    }
    
    
    ////////////////////////////////////////////////
    /////// LIBRARY OF CALCULATION FUNCTIONS ///////
    ////////////////////////////////////////////////
    
    double DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForce(SphericParticle* const element1,
                                                                    SphericParticle* const element2){
                
        //Get equivalent Radius
        const double my_radius     = element1->GetRadius();
        const double other_radius  = element2->GetRadius();
        const double cohesion_radius = my_radius > other_radius ? other_radius : my_radius;
        
        //Get cohesion Area
        const double cohesion_area = KRATOS_M_PI * cohesion_radius * cohesion_radius;
        
        //Get cohesion stress (user input)
        const double cohesion = 0.5 * ( element1->GetParticleCohesion() + element2->GetParticleCohesion() ); 
        
        //Force
        return cohesion_area * cohesion;
    }
    
    double DEMDiscontinuumConstitutiveLaw::CalculateStandardCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall){
        
        const double cohesion         = element->GetParticleCohesion();
        const double equiv_cohesion   = 0.5 * (cohesion + wall->GetProperties()[WALL_COHESION]);
        const double element_radius   = element->GetRadius();
        const double area             = KRATOS_M_PI * element_radius * element_radius;
        return equiv_cohesion * area;
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateStandardTangentialForce(const double normal_force,
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element1,
                                                    SphericParticle* const element2) {        
                        
        LocalElasticContactForce[0] += - mKt * LocalDeltDisp[0];  // 0: first tangential
        LocalElasticContactForce[1] += - mKt * LocalDeltDisp[1];  // 1: second tangential
        
        //Get equivalent tangent of friction angle
        const double my_tg_of_friction_angle        = element1->GetTgOfFrictionAngle();
        const double other_tg_of_friction_angle     = element2->GetTgOfFrictionAngle();
        const double equiv_tg_of_fri_ang            = 0.5 * (my_tg_of_friction_angle + other_tg_of_friction_angle);                
        double ShearForceNow = DEM_MODULUS_3(LocalElasticContactForce);
        double Frictional_ShearForceMax = equiv_tg_of_fri_ang * normal_force;

        if (Frictional_ShearForceMax < 0.0) {
            Frictional_ShearForceMax = 0.0;
        }

        if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
            double ShearForceNowMaxRatio = Frictional_ShearForceMax / ShearForceNow;
            LocalElasticContactForce[0] = ShearForceNowMaxRatio * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = ShearForceNowMaxRatio * LocalElasticContactForce[1];
            sliding = true;
        }        
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateStandardTangentialForceWithFEM(const double OldLocalContactForce[3],
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall) {        
        
        const double WallBallFriction = wall->mTgOfFrictionAngle;         //TODO: CHECK THIS (SHOULD IT BE AN AVERAGE?)      
        
        LocalElasticContactForce[0] = OldLocalContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalContactForce[1] - mKt * LocalDeltDisp[1];

        const double ShearForceMax = LocalElasticContactForce[2] * WallBallFriction;
        const double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        if (ShearForceNow > ShearForceMax) {
            double inv_ShearForceNow = 1.0 / ShearForceNow;
            LocalElasticContactForce[0] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[1];
            sliding = true;
        }                     
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateStandardViscoDampingForce(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element1,
                                                                SphericParticle* const element2){                                        
        
        //Get equivalent mass
        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
        const double equiv_mass = sqrt(my_mass * other_mass);
                
        //Get equivalent ln of restitution coefficient
        const double my_ln_of_restit_coeff    = element1->GetLnOfRestitCoeff();
        const double other_ln_of_restit_coeff = element2->GetLnOfRestitCoeff();
        const double equiv_ln_of_restit_coeff = 0.5 * (my_ln_of_restit_coeff + other_ln_of_restit_coeff);
        double equiv_visco_damp_coeff_normal  = 0.0;
        
        if (my_ln_of_restit_coeff > 0.0 || other_ln_of_restit_coeff > 0.0) { // Limit expressions when the restitution coefficient tends to 0. Variable lnRestitCoeff is set to 1.0 (instead of minus infinite) by the problem type.
            equiv_visco_damp_coeff_normal = 2.0 * sqrt(equiv_mass * mKn);
        }

        else {
            equiv_visco_damp_coeff_normal = - 2.0 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * mKn / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + KRATOS_M_PI * KRATOS_M_PI));
        }

        const double aux_norm_to_tang = sqrt(mKt / mKn);
        const double equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;                
        
        //double RCTS = 2.0 * sqrt(my_mass / mKn);
        //std::cout << "Rayleigh's critical time step is = " << RCTS << std::endl;
        //std::cout << "1% of Rayleigh's critical time step is = " << 0.01 * RCTS << std::endl;
        
        //DAMPING FORCE
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
        
        if ( !sliding ) { 
            ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }
    }
    
    void DEMDiscontinuumConstitutiveLaw::CalculateStandardViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                bool sliding,
                                                                SphericParticle* const element,
                                                                DEMWall* const wall){                                        
        
        double equiv_visco_damp_coeff_normal     = 0.0;
        const double my_mass    = element->GetMass();
        const double my_ln_of_restit_coeff = element->GetLnOfRestitCoeff();

        if (my_ln_of_restit_coeff > 0.0) { equiv_visco_damp_coeff_normal = 2.0 * sqrt(my_mass * mKn); }

        else {
            equiv_visco_damp_coeff_normal = -2.0 * my_ln_of_restit_coeff * sqrt(my_mass * mKn / (my_ln_of_restit_coeff * my_ln_of_restit_coeff + KRATOS_M_PI * KRATOS_M_PI));
        }

        const double aux_norm_to_tang = sqrt(mKt / mKn);
        const double equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;
        
        //double RCTS = 2.0 * sqrt(my_mass / mKn);
        //std::cout << "Rayleigh's critical time step is = " << RCTS << std::endl;
        //std::cout << "1% of Rayleigh's critical time step is = " << 0.01 * RCTS << std::endl;
                
        //DAMPING FORCE
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
        
        if ( !sliding ) { 
            ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }
    }
            
    void DEMDiscontinuumConstitutiveLaw::CalculateForces(double LocalElasticContactForce[3],
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

    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
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
    
    
    void DEMDiscontinuumConstitutiveLaw::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
            double &equiv_visco_damp_coeff_tangential,
            SphericParticle* element1,
            SphericParticle* element2,
            double kn_el,
            double kt_el) {
        
        KRATOS_TRY 
        double aux_norm_to_tang = 0.0;
        const double mRealMass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        const double mDempack_local_damping = element1->GetProperties()[DEMPACK_LOCAL_DAMPING];

        equiv_visco_damp_coeff_normal = mDempack_local_damping * 2.0 * sqrt(kn_el / (mRealMass + other_real_mass)) * (sqrt(mRealMass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // Dempack no ho fa servir...
        KRATOS_CATCH("")  
    }   

} // KRATOS
