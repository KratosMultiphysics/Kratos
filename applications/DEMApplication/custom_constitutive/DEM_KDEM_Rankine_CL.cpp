// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Rankine_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Rankine::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Rankine(*this));
        return p_clone;
    }

    void DEM_KDEM_Rankine::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Rankine to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_Rankine::Check(Properties::Pointer pProp) const {
        DEM_KDEM::Check(pProp);
        if(!pProp->Has(CONTACT_SIGMA_MIN)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CONTACT_SIGMA_MIN should be present in the properties when using DEM_KDEM_Rankine. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CONTACT_SIGMA_MIN) = 0.0;
        }
    }

    void DEM_KDEM_Rankine::CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2){

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (failure_type == 0) {
            double tension_limit = 0.5 * (GetContactSigmaMax(element1) + GetContactSigmaMax(element2)); //N/m2

            BoundedMatrix<double, 3, 3> average_stress_tensor = ZeroMatrix(3,3);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
                }
            }

            Vector principal_stresses(3);
            noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);

            for (int i=0; i<3; i++) {
                if(principal_stresses[i] > tension_limit) {
                    failure_type = 4;
                    break;
                }
            }
        }

    }

    void DEM_KDEM_Rankine::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        //Firstly, we check that the bond is not broken (it can break in any state of forces or indentations, because breakage depends on the stress tensor)
        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (indentation >= 0.0) { //COMPRESSION This response is the same for broken or intact bonds!
            LocalElasticContactForce[2] = kn_el * indentation;
        }
        else {
            if (failure_type > 0) {
                LocalElasticContactForce[2] = 0.0;
            }
            else {
                LocalElasticContactForce[2] = kn_el * indentation;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Rankine::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
            double& contact_tau,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            DenseVector<int>& search_control_vector,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential

        if (failure_type == 0) {
            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }
        }
        else {
            LocalElasticExtraContactForce[0] = 0.0;
            LocalElasticExtraContactForce[1] = 0.0;

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            const double equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());

            if(equiv_tg_of_fri_ang < 0.0) {
                KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
            }

            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }

    bool DEM_KDEM_Rankine::CheckRequirementsOfStressTensor() {

        return true;

    }

} // namespace Kratos
