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
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template<class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateMaterialResponsePK1(
        ConstitutiveLaw::Parameters& rValues
        )
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, r_deformation_gradient_f, determinant_f,
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

    Vector& r_stress_vector                = rValues.GetStressVector();
    const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(r_stress_vector, r_deformation_gradient_f, determinant_f,
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
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector         = rValues.GetStrainVector();

    // We get the Deformation gradient F
    const Matrix& rF = rValues.GetDeformationGradientF();
    const Matrix& r_B = prod(rF, trans(rF));
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(r_B, r_strain_vector);

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
                noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            } else {
                noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            }
        }
    } else { // We check for plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain = this->GetPlasticStrain();

            BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
            BoundedArrayType predictive_stress_vector;
            noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
            BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
            BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);

            // Compute the plastic parameters
            const double F = TConstLawIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

            if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
                noalias(r_integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                TConstLawIntegratorType::IntegrateStressVector(
                    predictive_stress_vector, r_strain_vector, uniaxial_stress,
                    threshold, plastic_denominator, f_flux, g_flux,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, plastic_strain, rValues,
                    characteristic_length);
                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_Kirchhoff, plastic_strain);
                }
            }
        }
    }
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

    Vector& r_stress_vector       = rValues.GetStressVector();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f  = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    r_stress_vector       /= determinant_f;
    r_constitutive_matrix /= determinant_f;
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
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector         = rValues.GetStrainVector();

    // We get the Deformation gradient F
    const Matrix& rF = rValues.GetDeformationGradientF();
    const Matrix& r_B = prod(rF, trans(rF));
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(r_B, r_strain_vector);

    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // We compute the stress or the constitutive matrix
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

        // We get some variables
        double& r_threshold = this->GetThreshold();
        double& r_plastic_dissipation = this->GetPlasticDissipation();
        Vector& r_plastic_strain = this->GetPlasticStrain();

        BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
        BoundedArrayType predictive_stress_vector;
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - r_plastic_strain);

        // Initialize Plastic Parameters
        double uniaxial_stress = 0.0, plastic_denominator = 0.0;
        BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
        BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);

        // Compute the plastic parameters
        const double F = TConstLawIntegratorType::CalculatePlasticParameters(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            r_threshold, plastic_denominator, f_flux, g_flux,
            r_plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, rValues, characteristic_length,
            r_plastic_strain);

        if (F > std::abs(1.0e-4 * r_threshold)) { // Plastic case
            // While loop backward euler
            /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                r_threshold, plastic_denominator, f_flux, g_flux,
                r_plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, r_plastic_strain, rValues,
                characteristic_length);
        }
    }
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
void GenericFiniteStrainIsotropicPlasticity<TConstLawIntegratorType>::
    CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure& rStressMeasure,
        const Vector& rPlasticStrain
        )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    // if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
    //     // Already stored in rValues.GetConstitutiveMatrix()...
    // } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
    //     // Calculates the Tangent Constitutive Tensor by perturbation (first order)
    //     TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 1);
    // } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
    //     // Calculates the Tangent Constitutive Tensor by perturbation (second order)
    //     TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 2);
    // } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {
    //     const Vector num = prod(rValues.GetConstitutiveMatrix(), rPlasticStrain);
    //     const double denom = inner_prod(rValues.GetStrainVector(), num);
    //     noalias(rValues.GetConstitutiveMatrix()) -= outer_prod(num, num) / denom;
    // }
    BaseType::CalculateElasticMatrix( rValues.GetConstitutiveMatrix(), rValues);
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
    if (rThisVariable == UNIAXIAL_STRESS) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);
        const Vector& r_stress_vector = rParameterValues.GetStressVector();
        const Vector& r_strain_vector = rParameterValues.GetStrainVector();

        BoundedArrayType aux_stress_vector = r_stress_vector;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(aux_stress_vector,
            r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
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
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == HENCKY_STRAIN_VECTOR ||
        rThisVariable == BIOT_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the strain
        if (rThisVariable == STRAIN) {
            // HyperElasticIsotropicKirchhoff3D::CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
                Matrix identity_matrix(Dimension, Dimension);
                noalias(identity_matrix) = IdentityMatrix(Dimension);
                const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
                const Matrix E_matrix = 0.5 * (prod(trans(deformation_gradient_f), deformation_gradient_f) - identity_matrix);
                Vector& r_strain_vector = rParameterValues.GetStrainVector();
                noalias(r_strain_vector) = MathUtils<double>::StrainTensorToVector(E_matrix, VoigtSize);
        } else if (rThisVariable == ALMANSI_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix B_tensor = prod(deformation_gradient_f, trans(deformation_gradient_f));
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateAlmansiStrain(B_tensor, r_strain_vector);
        } else if (rThisVariable == HENCKY_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_tensor, r_strain_vector);
        } else if (rThisVariable == BIOT_STRAIN_VECTOR) {
            const Matrix& deformation_gradient_f = rParameterValues.GetDeformationGradientF();
            const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateBiotStrain(C_tensor, r_strain_vector);
        }

        rValue = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        if (rThisVariable == STRESSES) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        } if (rThisVariable == CAUCHY_STRESS_VECTOR) {
            this->CalculateMaterialResponseCauchy(rParameterValues);
        } if (rThisVariable == PK2_STRESS_VECTOR) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        }
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
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
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(this->GetPlasticStrain());
    } else if (this->Has(rThisVariable)) {
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
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
        ) const
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
