// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//  Collaborator:
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_advanced_constitutive/generic_small_strain_kinematic_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
                noalias(r_stress_vector) = prod( r_constitutive_matrix, r_strain_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

        //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
        if ( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateValue(rValues, STRAIN, r_strain_vector);
        }

        // Elastic Matrix
        if ( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
        }

        // We compute the stress
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            // Elastic Matrix
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain      = this->GetPlasticStrain();
            Vector back_stress_vector  = this->GetBackStressVector();
            const Vector previous_stress_vector = this->GetPreviousStressVector();

            array_1d<double, VoigtSize> predictive_stress_vector, kin_hard_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::U_P_LAW)) {
                predictive_stress_vector = rValues.GetStressVector();
            } else {
                // S0 = r_constitutive_matrix:(E-Ep)
                predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
            }

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
            array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

            // Kinematic back stress substracted
            noalias(kin_hard_stress_vector) = predictive_stress_vector - back_stress_vector;

            const double F = TConstLawIntegratorType::CalculatePlasticParameters(
                kin_hard_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain, back_stress_vector);

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
                    characteristic_length, back_stress_vector,
                    previous_stress_vector);

                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
                }
            }
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if ( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // We get some variables
    double threshold = this->GetThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain      = this->GetPlasticStrain();
    Vector back_stress_vector  = this->GetBackStressVector();
    const Vector previous_stress_vector = this->GetPreviousStressVector();

    array_1d<double, VoigtSize> predictive_stress_vector, kin_hard_stress_vector;
    if (r_constitutive_law_options.Is(ConstitutiveLaw::U_P_LAW)) {
        predictive_stress_vector = rValues.GetStressVector();
    } else {
        // S0 = r_constitutive_matrix:(E-Ep)
        predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
    }

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
    array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
    array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

    // Kinematic back stress substracted
    noalias(kin_hard_stress_vector) = predictive_stress_vector - back_stress_vector;

    const double threshold_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
        kin_hard_stress_vector, r_strain_vector, uniaxial_stress,
        threshold, plastic_denominator, f_flux, g_flux,
        plastic_dissipation, plastic_strain_increment,
        r_constitutive_matrix, rValues, characteristic_length,
        plastic_strain, back_stress_vector);

    if (threshold_indicator > std::abs(1.0e-4 * threshold)) {
        // while loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, plastic_strain, rValues,
            characteristic_length, back_stress_vector,
            previous_stress_vector);
    }

    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

    mPlasticDissipation = plastic_dissipation;
    mThreshold = threshold;

    noalias(mPlasticStrain) = plastic_strain;
    noalias(mPreviousStressVector) = predictive_stress_vector;
    noalias(mBackStressVector) = back_stress_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        return true;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        return true;
    } else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mPlasticDissipation = rValue[0];
        for (std::size_t i=0; i < VoigtSize; ++i)
            mPlasticStrain[i] = rValue[i + 1];
        return;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
        return;
    }
    BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        rValue.resize(1 + VoigtSize);
        rValue[0] = mPlasticDissipation;
        for (std::size_t i=0; i < VoigtSize; ++i)
            rValue[i + 1] = mPlasticStrain[i];
        return rValue;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
        return rValue;
    }
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        rValue = MathUtils<double>::StressVectorToTensor(mPreviousStressVector);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
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
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( aux_stress_vector, r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
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

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );

        // Compute the equivalent plastic strain
        double uniaxial_stress;
        this->CalculateValue(rParameterValues, UNIAXIAL_STRESS, uniaxial_stress);
        TConstLawIntegratorType::CalculateEquivalentPlasticStrain(r_stress_vector, uniaxial_stress, mPlasticStrain, 0.0, rParameterValues, rValue);
        return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == BACK_STRESS_VECTOR) {
        const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rParameterValues.GetElementGeometry());
        const Flags& r_constitutive_law_options = rParameterValues.GetOptions();

        // We get the strain vector
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rParameterValues.GetConstitutiveMatrix();
        this->CalculateValue(rParameterValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
        if ( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateValue(rParameterValues, STRAIN, r_strain_vector);
        }

        // We get some variables
        double threshold = this->GetThreshold();
        double plastic_dissipation = this->GetPlasticDissipation();
        Vector plastic_strain      = this->GetPlasticStrain();
        Vector back_stress_vector  = this->GetBackStressVector();
        const Vector previous_stress_vector = this->GetPreviousStressVector();

        array_1d<double, VoigtSize> predictive_stress_vector, kin_hard_stress_vector;
        if (r_constitutive_law_options.Is(ConstitutiveLaw::U_P_LAW)) {
            predictive_stress_vector = rParameterValues.GetStressVector();
        } else {
            // S0 = r_constitutive_matrix:(E-Ep)
            predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        }

        // Initialize Plastic Parameters
        double uniaxial_stress = 0.0, plastic_denominator = 0.0;
        array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
        array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
        array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

        // Kinematic back stress substracted
        noalias(kin_hard_stress_vector) = predictive_stress_vector - back_stress_vector;

        const double threshold_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
            kin_hard_stress_vector, r_strain_vector, uniaxial_stress,
            threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, rParameterValues, characteristic_length,
            plastic_strain, back_stress_vector);

        if (threshold_indicator > std::abs(1.0e-4 * threshold)) {
            // while loop backward euler
            /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, plastic_strain, rParameterValues,
                characteristic_length, back_stress_vector,
                previous_stress_vector);
        }
        rValue = back_stress_vector;
    } else {
        BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == BACK_STRESS_TENSOR) {
        Vector aux_value(VoigtSize);
        this->CalculateValue(rParameterValues, BACK_STRESS_VECTOR, aux_value);
        rValue = MathUtils<double>::StressVectorToTensor(aux_value);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Check(
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

template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

} // namespace Kratos
