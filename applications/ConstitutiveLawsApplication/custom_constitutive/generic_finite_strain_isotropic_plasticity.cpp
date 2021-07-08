// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateMaterialResponsePK1(
        ConstitutiveLaw::Parameters& rValues
        )
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f,
        ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateMaterialResponsePK2(
        ConstitutiveLaw::Parameters& rValues
        )
{
    this->CalculateMaterialResponseKirchhoff(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f,
        ConstitutiveLaw::StressMeasure_Kirchhoff, ConstitutiveLaw::StressMeasure_PK2);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateMaterialResponseKirchhoff(
        ConstitutiveLaw::Parameters& rValues
        )
{
    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateMaterialResponseCauchy(
        ConstitutiveLaw::Parameters& rValues
        )
{
    this->CalculateMaterialResponseKirchhoff(rValues);

    Vector& stress_vector       = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f  = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    stress_vector       /= determinant_f;
    constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    FinalizeMaterialResponsePK1(
        ConstitutiveLaw::Parameters& rValues
        )
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    FinalizeMaterialResponsePK2(
        ConstitutiveLaw::Parameters& rValues
        )
{
    FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    FinalizeMaterialResponseKirchhoff(
        ConstitutiveLaw::Parameters& rValues
        )
{
    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    FinalizeMaterialResponseCauchy(
        ConstitutiveLaw::Parameters& rValues
        )
{
    FinalizeMaterialResponseKirchhoff(rValues);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
double& GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        )
{
    // if (rThisVariable == UNIAXIAL_STRESS) {
    //     // Get Values to compute the constitutive law:
    //     Flags& r_flags = rParameterValues.GetOptions();

    //     // Previous flags saved
    //     const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    //     const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

    //     r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

    //     // Calculate the stress vector
    //     CalculateMaterialResponseCauchy(rParameterValues);
    //     const Vector& r_stress_vector = rParameterValues.GetStressVector();
    //     const Vector& r_strain_vector = rParameterValues.GetStrainVector();

    //     BoundedArrayType aux_stress_vector = r_stress_vector;
    //     TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(aux_stress_vector,
    //         r_strain_vector, rValue, rParameterValues);

    //     // Previous flags restored
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    // } else if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
    //     // Calculate the stress vector
    //     Vector stress_vector;
    //     this->CalculateValue(rParameterValues, PK2_STRESS_VECTOR, stress_vector);

    //     // The plastic deformation gradient
    //     const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

    //     // We backup the deformation gradient
    //     const double& deformation_gradient_determinant_backup = rParameterValues.GetDeterminantF();
    //     const Matrix& r_deformation_gradient_backup = rParameterValues.GetDeformationGradientF();

    //     rParameterValues.SetDeterminantF(MathUtils<double>::Det(r_plastic_deformation_gradient));
    //     rParameterValues.SetDeformationGradientF(r_plastic_deformation_gradient);
    //     Vector plastic_strain;
    //     this->CalculateValue(rParameterValues, GREEN_LAGRANGE_STRAIN_VECTOR, plastic_strain);

    //     // Compute the equivalent plastic strain
    //     double uniaxial_stress;
    //     this->CalculateValue(rParameterValues, UNIAXIAL_STRESS, uniaxial_stress);
    //     TConstLawIntegratorType::CalculateEquivalentPlasticStrain(stress_vector,
    //         uniaxial_stress, plastic_strain, 0.0, rParameterValues, rValue);

    //     // We revert the deformation gradient
    //     rParameterValues.SetDeterminantF(deformation_gradient_determinant_backup);
    //     rParameterValues.SetDeformationGradientF(r_deformation_gradient_backup);

    //     return rValue;
    // } else {
    //     return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    // }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
Vector& GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
Matrix& GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
    )
{
    // if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
    //     const Matrix C_plastic_tensor = prod(trans( mPlasticDeformationGradient), mPlasticDeformationGradient);
    //     Vector plastic_strain(VoigtSize);
    //     AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_plastic_tensor, plastic_strain);
    //     rValue = MathUtils<double>::StrainVectorToTensor(plastic_strain);
    // } else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
    //     // The plastic deformation gradient
    //     const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

    //     // We backup the deformation gradient
    //     const Matrix& r_deformation_gradient_backup = rParameterValues.GetDeformationGradientF();

    //     // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
    //     Matrix inverse_F_p ( Dimension, Dimension );
    //     double aux_det_Fp = 0;
    //     MathUtils<double>::InvertMatrix( r_plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
    //     const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

    //     rParameterValues.SetDeformationGradientF(elastic_deformation_gradient);
    //     Vector elastic_stress;
    //     this->CalculateValue(rParameterValues, CAUCHY_STRESS_VECTOR, elastic_stress);

    //     rValue = MathUtils<double>::StressVectorToTensor(elastic_stress);

    //     // We revert the deformation gradient
    //     rParameterValues.SetDeformationGradientF(r_deformation_gradient_backup);

    //     return rValue;
    // } else if (this->Has(rThisVariable)) {
    //     return this->GetValue(rThisVariable, rValue);
    // } else {
    //     return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    // }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TConstLawIntegratorType>
int GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    const int check_integrator = TConstLawIntegratorType::Check(rMaterialProperties);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    if ((check_base + check_integrator) > 0) return 1;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

} // namespace Kratos
