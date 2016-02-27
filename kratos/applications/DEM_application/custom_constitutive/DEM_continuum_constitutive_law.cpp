#include "DEM_application.h"
#include "DEM_continuum_constitutive_law.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw() {
        //std::cout << " DEMContinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMContinuumConstitutiveLaw

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw) {
        //std::cout << " DEMContinuumConstitutiveLaw copy constructor.." << std::endl;
    }
    
    std::string DEMContinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMContinuumConstitutiveLaw::Initialize(const ProcessInfo& r_process_info) {
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEMContinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw() { //std::cout << "Law destructor..." ;
    }

    void DEMContinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
                                                            double ViscoDampingLocalContactForce[3],
                                                            double indentation,
                                                            double equiv_visco_damp_coeff_normal,
                                                            double equiv_visco_damp_coeff_tangential,
                                                            bool sliding,
                                                            int failure_id) {

        //*** component-wise since localContactForce and RelVel have in principle no relationship.
        // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
        // But in opposite direction the visco damping can't overpass the force...

        KRATOS_TRY  

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];

        if (sliding == false) { //only applied when no sliding to help the regularized friction law or the spring convergence
            ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        }

        KRATOS_CATCH("")
    }
    
    void DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                                       SphericContinuumParticle* neighbor,
                                                                       double equiv_young,
                                                                       double distance,
                                                                       double calculation_area,
                                                                       double LocalCoordSystem[3][3],
                                                                       array_1d<double, 3>& mContactMoment) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments) should not be called.","")
    }
    
    void DEMContinuumConstitutiveLaw::AddPoissonContribution(const double equiv_poisson, 
                                                            double LocalCoordSystem[3][3], 
                                                            double& normal_force, 
                                                            double calculation_area, Matrix* mSymmStressTensor){    
    }
    
} //kratos
