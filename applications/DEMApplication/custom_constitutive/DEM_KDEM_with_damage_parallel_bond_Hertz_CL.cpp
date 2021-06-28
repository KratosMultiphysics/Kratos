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

    void DEM_KDEM_with_damage_parallel_bond_Hertz::SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond_Hertz to Properties " << pProp->Id() <<" with given parameters"<< std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());

        TransferParametersToProperties(parameters, pProp);

        this->Check(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                             double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

        KRATOS_TRY

        // Unbonded elastic constants
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        double sqrt_equiv_radius_and_indentation;
        if (indentation > 0.0) {
            sqrt_equiv_radius_and_indentation = sqrt(equiv_radius * indentation);
        } else {
            sqrt_equiv_radius_and_indentation = 0.0;
        }

        //Normal and Tangent elastic constants
        mUnbondedNormalElasticConstant = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mUnbondedTangentialElasticConstant = 4.0 * equiv_shear * mUnbondedNormalElasticConstant / equiv_young;

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();
        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);
        const double& equiv_gamma   = (*mpProperties)[DAMPING_GAMMA];

        mUnbondedEquivViscoDampCoeffNormal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedNormalElasticConstant);
        mUnbondedEquivViscoDampCoeffTangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedTangentialElasticConstant);

        const double bonded_equiv_young = (*mpProperties)[BONDED_MATERIAL_YOUNG_MODULUS];
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

    double DEM_KDEM_with_damage_parallel_bond_Hertz::LocalMaxSearchDistance(const int i,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2) {

        double tension_limit;

        // calculation of equivalent Young modulus
        const double bonded_equivalent_young = (*mpProperties)[BONDED_MATERIAL_YOUNG_MODULUS];

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0.0;

        Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);
        GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);

        double radius_sum = my_radius + other_radius;
        double initial_delta = element1->GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants
        double kn_el = bonded_equivalent_young * calculation_area / initial_dist;

        tension_limit = GetContactSigmaMax();

        const double Ntstr_el = tension_limit * calculation_area;
        double u1 = Ntstr_el / kn_el;
        if (u1 > 2.0 * radius_sum) {
            u1 = 2.0 * radius_sum;
        } // avoid error in special cases with too high tensile
        return u1;
    }
} // namespace Kratos
