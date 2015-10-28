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

    void DEMContinuumConstitutiveLaw::Initialize(const ProcessInfo& rCurrentProcessInfo) {
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

} //kratos
