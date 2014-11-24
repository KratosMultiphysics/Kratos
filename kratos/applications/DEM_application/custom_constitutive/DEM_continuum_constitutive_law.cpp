#include "DEM_application.h"
#include "DEM_continuum_constitutive_law.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw() {
        //std::cout << " DEMContinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMContinuumConstitutiveLaw

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw) {
        //std::cout << " DEMContinuumConstitutiveLaw copy constructor.." << std::endl;
    }

    void DEMContinuumConstitutiveLaw::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << " Assigning DEMContinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw() { //std::cout << "Law destructor..." ;
    }

    void DEMContinuumConstitutiveLaw::CalculateContactForces(double mRadius,
            double mSqrtOfRealMass,
            double other_radius,
            double otherSqrtMass,
            double distance,
            double initial_delta,
            int& neighbour_failure_id,
            ProcessInfo& rCurrentProcessInfo,
            PropertiesProxy *myProperties,
            PropertiesProxy *neighbourProperties,
            int mapping_new_ini,
            int mapping_new_cont,
            unsigned int i_neighbour_count,
            double LocalElasticContactForce[3],
            double ViscoDampingLocalContactForce[3],
            double LocalDeltDisp[3],
            Vector mcont_ini_neigh_area,
            array_1d<double, 4 > &mHistory_mapping_new_cont,
            double mDempack_damping,
            int mDampType,
            int mIniNeighbourFailureId_mapping_new_ini,
            double LocalCoordSystem[3][3],
            double RelVel[3]) {
    }

    void DEMContinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
            double ViscoDampingLocalContactForce[3],
            double indentation,
            double equiv_visco_damp_coeff_normal,
            double equiv_visco_damp_coeff_tangential,
            bool sliding,
            int mDampType) {

        // The comprovation is component-wise since localContactForce and RelVel have in principle no relationship.
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
    }    

} //kratos
