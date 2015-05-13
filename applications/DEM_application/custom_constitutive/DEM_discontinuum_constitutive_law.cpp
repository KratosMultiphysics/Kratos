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
