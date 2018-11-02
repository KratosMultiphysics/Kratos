// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress plasticity
    Vector& integrated_stress_vector = rValues.GetStressVector();
    Matrix& tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const double characteristic_length = rValues.GetElementGeometry().Length();
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    // We compute the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

        // We get some variables
<<<<<<< HEAD
        double threshold = this->GetThreshold();
        double plastic_dissipation = this->GetPlasticDissipation();
        Vector plastic_strain = this->GetPlasticStrain();
=======
        double r_threshold = this->GetThreshold();
        double r_plastic_dissipation = this->GetPlasticDissipation();
        Vector r_plastic_strain = this->GetPlasticStrain();
>>>>>>> master

        array_1d<double, VoigtSize> predictive_stress_vector;
        if( r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW ) ) {
            predictive_stress_vector = rValues.GetStressVector();
        } else {
            // S0 = r_constitutive_matrix:(E-Ep)
            predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        }

        // Initialize Plastic Parameters
        double uniaxial_stress = 0.0, plastic_denominator = 0.0;
        array_1d<double, VoigtSize> f_flux(VoigtSize, 0.0); // DF/DS
        array_1d<double, VoigtSize> g_flux(VoigtSize, 0.0); // DG/DS
        array_1d<double, VoigtSize> plastic_strain_increment(VoigtSize, 0.0);

        TConstLawIntegratorType::CalculatePlasticParameters(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, rValues, characteristic_length,
            plastic_strain);

        const double F = uniaxial_stress - threshold;

        if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
            noalias(integrated_stress_vector) = predictive_stress_vector;
<<<<<<< HEAD
            this->SetNonConvPlasticDissipation(plastic_dissipation);
            this->SetNonConvPlasticStrain(plastic_strain);
            this->SetNonConvThreshold(threshold);
=======
            this->SetNonConvPlasticDissipation(r_plastic_dissipation);
            this->SetNonConvPlasticStrain(r_plastic_strain);
            this->SetNonConvThreshold(r_threshold);
>>>>>>> master
            this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(tangent_tensor) = r_constitutive_matrix;
            }
        } else { // Plastic case
            // while loop backward euler
            /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
            TConstLawIntegratorType::IntegrateStressVector(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, plastic_strain, rValues,
                characteristic_length);
            noalias(integrated_stress_vector) = predictive_stress_vector;

<<<<<<< HEAD
            this->SetNonConvPlasticDissipation(plastic_dissipation);
            this->SetNonConvPlasticStrain(plastic_strain);
            this->SetNonConvThreshold(threshold);
=======
            this->SetNonConvPlasticDissipation(r_plastic_dissipation);
            this->SetNonConvPlasticStrain(r_plastic_strain);
            this->SetNonConvThreshold(r_threshold);
>>>>>>> master
            this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
                noalias(tangent_tensor) = rValues.GetConstitutiveMatrix();
            }
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::InitializeMaterial(
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
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->SetPlasticDissipation(this->GetNonConvPlasticDissipation());
    this->SetThreshold(this->GetNonConvThreshold());
    this->SetPlasticStrain(this->GetNonConvPlasticStrain());
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        mUniaxialStress = rValue;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        rValue = mUniaxialStress;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
        const Vector& r_stress_vector = rParameterValues.GetStressVector();
        TConstLawIntegratorType::CalculateEquivalentPlasticStrain(r_stress_vector, mUniaxialStress, mPlasticStrain, 0.0, rParameterValues, rValue);

        return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        //1.-Compute total deformation gradient
        const Matrix& DeformationGradientF = rParameterValues.GetDeformationGradientF();
        //2.-Right Cauchy-Green tensor C
        const Matrix right_cauchy_green = prod(trans(DeformationGradientF), DeformationGradientF);
        Vector strain_vector = ZeroVector(6);

        //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
        strain_vector[0] = 0.5 * (right_cauchy_green(0, 0) - 1.00);
        strain_vector[1] = 0.5 * (right_cauchy_green(1, 1) - 1.00);
        strain_vector[2] = 0.5 * (right_cauchy_green(2, 2) - 1.00);
        strain_vector[3] = right_cauchy_green(0, 1); // xy
        strain_vector[4] = right_cauchy_green(1, 2); // yz
        strain_vector[5] = right_cauchy_green(0, 2); // xz

        Matrix constitutive_matrix;
        this->CalculateElasticMatrix(constitutive_matrix, rParameterValues);

        array_1d<double,VoigtSize> tmp = prod(constitutive_matrix, strain_vector - mPlasticStrain);
        rValue = MathUtils<double>::StressVectorToTensor(tmp);
        return rValue;
    } else if (this->Has(rThisVariable)) {
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>::Check(
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

template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;

} // namespace Kratos
