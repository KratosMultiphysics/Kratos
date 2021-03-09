// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_parallel_bond_Hertz_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage_parallel_bond_Hertz::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage_parallel_bond_Hertz(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond_Hertz to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                             double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

        KRATOS_TRY

        //TODO: Sometimes we do not compute mean values this way. Sometimes we use 2xy/(x+y)
        const double unbonded_equivalent_young = 0.5 * (element1->GetProperties()[LOOSE_MATERIAL_YOUNG_MODULUS] + element2->GetProperties()[LOOSE_MATERIAL_YOUNG_MODULUS]);
        const double unbonded_equivalent_shear = unbonded_equivalent_young / (2.0 * (1 + equiv_poisson));

        // Unbonded elastic constants
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_radius * indentation);
        mUnbondedNormalElasticConstant = 2.0 * unbonded_equivalent_young * sqrt_equiv_radius_and_indentation;
        mUnbondedTangentialElasticConstant = 4.0 * unbonded_equivalent_shear * mUnbondedNormalElasticConstant / equiv_young;

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        mEquivViscoDampCoeffNormal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedNormalElasticConstant);
        mEquivViscoDampCoeffTangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedTangentialElasticConstant);
        //

        const double bonded_equiv_young = equiv_young - unbonded_equivalent_young;
        const double bonded_equiv_shear = bonded_equiv_young / (2.0 * (1 + equiv_poisson));

        kn_el = bonded_equiv_young * calculation_area / initial_dist;
        kt_el = bonded_equiv_shear * calculation_area / initial_dist;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz::ComputeNormalUnbondedForce(double indentation) {

        KRATOS_TRY

        mUnbondedLocalElasticContactForce2 = 0.666666666666666666667 * mUnbondedNormalElasticConstant * indentation;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz::CalculateViscoDamping(double LocalRelVel[3],
                                         double ViscoDampingLocalContactForce[3],
                                         double indentation,
                                         double equiv_visco_damp_coeff_normal,
                                         double equiv_visco_damp_coeff_tangential,
                                         bool& sliding,
                                         int failure_id) {

        KRATOS_TRY

        if (indentation > 0 && !sliding) {
            ViscoDampingLocalContactForce[0] = -mEquivViscoDampCoeffTangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] = -mEquivViscoDampCoeffTangential * LocalRelVel[1];
            ViscoDampingLocalContactForce[2] = -mEquivViscoDampCoeffNormal * LocalRelVel[2];
        }

        mViscoDampingLocalContactForce[0] = ViscoDampingLocalContactForce[0];
        mViscoDampingLocalContactForce[1] = ViscoDampingLocalContactForce[1];
        mViscoDampingLocalContactForce[2] = ViscoDampingLocalContactForce[2];


        if (!failure_id) {
            ViscoDampingLocalContactForce[0] -= equiv_visco_damp_coeff_tangential * LocalRelVel[0];
            ViscoDampingLocalContactForce[1] -= equiv_visco_damp_coeff_tangential * LocalRelVel[1];
            ViscoDampingLocalContactForce[2] -= equiv_visco_damp_coeff_normal * LocalRelVel[2];
        }

        KRATOS_CATCH("")
    }
} // namespace Kratos
