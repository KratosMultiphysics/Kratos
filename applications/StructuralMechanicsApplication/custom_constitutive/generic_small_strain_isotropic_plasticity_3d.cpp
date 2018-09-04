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
#include "custom_constitutive/generic_small_strain_isotropic_plasticity_3d.h"
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
template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
    // Integrate Stress plasticity
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const int voigt_size = this->GetStrainSize();
    Vector& integrated_stress_vector = rValues.GetStressVector();
    Matrix& tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const double characteristic_length = rValues.GetElementGeometry().Length();
    const Flags& r_constitututive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitututive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }

    // Elastic Matrix
    if( r_constitututive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    // We compute the stress
    if( r_constitututive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

        // We get some variables
        double& r_threshold = this->GetThreshold();
        double& r_plastic_dissipation = this->GetPlasticDissipation();
        Vector& r_plastic_strain = this->GetPlasticStrain();

        // S0 = r_constitutive_matrix:(E-Ep)
        Vector predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - r_plastic_strain);

        // Initialize Plastic Parameters
        double uniaxial_stress = 0.0, plastic_denominator = 0.0;
        Vector f_flux = ZeroVector(voigt_size), g_flux = ZeroVector(voigt_size); // DF/DS & DG/DS
        Vector plastic_strain_increment = ZeroVector(voigt_size);

        ConstLawIntegratorType::CalculatePlasticParameters(predictive_stress_vector, r_strain_vector,
                                                        uniaxial_stress, r_threshold, plastic_denominator, f_flux, g_flux, r_plastic_dissipation,
                                                        plastic_strain_increment, r_constitutive_matrix, r_material_properties, characteristic_length);

        const double F = uniaxial_stress - r_threshold;

        if (F <= std::abs(1.0e-4 * r_threshold)) { // Elastic case
            noalias(integrated_stress_vector) = predictive_stress_vector;
            this->SetNonConvPlasticDissipation(r_plastic_dissipation);
            this->SetNonConvPlasticStrain(r_plastic_strain);
            this->SetNonConvThreshold(r_threshold);

            if (r_constitututive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(tangent_tensor) = r_constitutive_matrix;
                this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
            }
        } else { // Plastic case
            // while loop backward euler
            /* Inside "IntegrateStressVector" the predictive_stress_vector
                is updated to verify the yield criterion */
            ConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, r_strain_vector,
                                                        uniaxial_stress, r_threshold, plastic_denominator, f_flux, g_flux,
                                                        r_plastic_dissipation, plastic_strain_increment, r_constitutive_matrix, r_plastic_strain,
                                                        r_material_properties, characteristic_length);
            noalias(integrated_stress_vector) = predictive_stress_vector;

            this->SetNonConvPlasticDissipation(r_plastic_dissipation);
            this->SetNonConvPlasticStrain(r_plastic_strain);
            this->SetNonConvThreshold(r_threshold);

            if (r_constitututive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
                noalias(tangent_tensor) = rValues.GetConstitutiveMatrix();
                this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
            }
        }
        // this->SetValue(INTEGRATED_STRESS_TENSOR, MathUtils<double>::StressVectorToTensor(integrated_stress_vector), rValues.GetProcessInfo());
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues)
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& r_material_properties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    this->SetPlasticDissipation(this->GetNonConvPlasticDissipation());
    this->SetThreshold(this->GetNonConvThreshold());
    this->SetPlasticStrain(this->GetNonConvPlasticStrain());
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponse(
    ConstitutiveLaw::Parameters &rValues,
    const StressMeasure &rStressMeasure)
{
} // FinalizeMaterialResponse

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::Has(const Variable<double> &rThisVariable)
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
bool GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::Has(const Variable<Vector> &rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::SetValue(
    const Variable<double> &rThisVariable,
    const double& rValue,
    const ProcessInfo &rCurrentProcessInfo)
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

template <class ConstLawIntegratorType>
void GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::SetValue(
    const Variable<Vector> &rThisVariable,
    const Vector& rValue,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double& GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<double> &rThisVariable,
    double& rValue)
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

template <class ConstLawIntegratorType>
Vector& GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<Vector> &rThisVariable,
    Vector& rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::GetValue(
    const Variable<Matrix> &rThisVariable,
    Matrix& rValue)
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
double& GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<double> &rThisVariable,
    double& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters &rParameterValues,
    const Variable<Matrix> &rThisVariable,
    Matrix& rValue)
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

        Matrix C;
        this->CalculateElasticMatrix(C, rParameterValues);

        rValue = MathUtils<double>::StressVectorToTensor(prod(C, strain_vector - mPlasticStrain));
        return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class ConstLawIntegratorType>
int GenericSmallStrainIsotropicPlasticity3D<ConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const int check_integrator = ConstLawIntegratorType::Check(rMaterialProperties);

    if ((check_base + check_integrator) > 0) return 1;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>;
template class GenericSmallStrainIsotropicPlasticity3D<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential>>>;

} // namespace Kratos
