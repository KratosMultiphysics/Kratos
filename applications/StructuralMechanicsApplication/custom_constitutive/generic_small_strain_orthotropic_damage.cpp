// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/generic_small_strain_orthotropic_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

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
template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Integrate Stress Damage
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Internal variables initialization
        Vector damages = mDamages;
        Vector thresholds = mThresholds;

        // S0 = C:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);

        // Now we compute the principal stresses
        array_1d<double, Dimension> principal_stresses_vector;
        ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stresses_vector, predictive_stress_vector);

        double uniaxial_stress = 0.0;
        // Now we compute the damages on each direction...
        for (unsigned int i = 0; i < Dimension; i++) {
            if (principal_stresses_vector[i] > tolerance) { // Damage only if positive
                TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
            }
            const double F = uniaxial_stress - thresholds[i];
            if (F > tolerance) { // Damage Case
                const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
                TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, damages[i], thresholds[i], rValues, characteristic_length);
            } else { // Elastic Case
                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    // noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;
                }
            }

        }



    }



} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateTangentTensor(
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
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::InitializeMaterial(
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

    Vector initial_thresholds(Dimension);
    noalias(initial_thresholds) = ZeroVector(Dimension);
    for (unsigned int i = 0; i < Dimension; i++)
        initial_thresholds[i] = initial_threshold;

    this->SetThresholds(initial_thresholds);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{






}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::Has(
    const Variable<double>& rThisVariable
    )
{
    // if (rThisVariable == DAMAGE) {
    //     return true;
    // } else if (rThisVariable == THRESHOLD) {
    //     return true;
    // } else if (rThisVariable == UNIAXIAL_STRESS) {
    //     return true;
    // } else {
    //     return BaseType::Has(rThisVariable);
    // }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::Has(
    const Variable<Vector>& rThisVariable
    )
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::Has(
    const Variable<Matrix>& rThisVariable
    )
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // if (rThisVariable == DAMAGE) {
    //     mDamage = rValue;
    // } else if (rThisVariable == THRESHOLD) {
    //     mThreshold = rValue;
    // } else if (rThisVariable == UNIAXIAL_STRESS) {
    //     mUniaxialStress = rValue;
    // } else {
    //     return BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // if (rThisVariable == DAMAGE) {
    //     rValue = mDamage;
    // } else if (rThisVariable == THRESHOLD) {
    //     rValue = mThreshold;
    // } else if (rThisVariable == UNIAXIAL_STRESS) {
    //     rValue = mUniaxialStress;
    // } else {
    //     return BaseType::GetValue(rThisVariable, rValue);
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateValue(
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
Matrix& GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // if (this->Has(rThisVariable)) {
    //     return this->GetValue(rThisVariable, rValue);
    // } else {
    //     return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    // }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::Check(
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

template <class TConstLawIntegratorType>
void GenericSmallStrainOrthotropicDamage<TConstLawIntegratorType>::CalculateSecantTensor(
    Matrix& rSecantTensor, 
    ConstitutiveLaw::Parameters& rValues, 
    const Vector& rDamages
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double nu = r_material_properties[POISSON_RATIO];

    if (rSecantTensor.size1() == VoigtSize)
        rSecantTensor.resize(VoigtSize, VoigtSize);
    noalias(rSecantTensor) = ZeroMatrix(VoigtSize, VoigtSize);

    if (Dimension == 3) { // 3D version
        const double c1 = E / ((1.00 + nu) * (1 - 2 * nu));
        const double c2 = c1 * (1 - nu); // Cii
        const double c3 = c1 * nu;  // Cij
        const double c4 = c1 * 0.5 * (1 - 2 * nu); // Gij
        rSecantTensor(0,0) = (1.0 - rDamages[0]) * c2;
        rSecantTensor(1,1) = (1.0 - rDamages[1]) * c2;
        rSecantTensor(2,2) = (1.0 - rDamages[2]) * c2;
        rSecantTensor(0,1) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c3;
        rSecantTensor(0,2) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[2]))*c3;
        rSecantTensor(1,0) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c3;
        rSecantTensor(1,2) = std::sqrt((1.0 - rDamages[1])*(1.0 - rDamages[2]))*c3;
        rSecantTensor(2,0) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[2]))*c3;
        rSecantTensor(2,1) = std::sqrt((1.0 - rDamages[1])*(1.0 - rDamages[2]))*c3;
        rSecantTensor(3,3) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c4;
        rSecantTensor(4,4) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[2]))*c4;
        rSecantTensor(5,5) = std::sqrt((1.0 - rDamages[1])*(1.0 - rDamages[2]))*c4;
    } else { // 2D Case -> Plane strain
        const double c0 = E / ((1.00 + nu)*(1 - 2 * nu));
        const double c1 = (1.00 - nu)*c0;
        const double c2 = c0 * nu;
        const double c3 = (0.5 - nu)*c0;

        rSecantTensor(0,0) = (1.0 - rDamages[0])*c1;
        rSecantTensor(1,1) = (1.0 - rDamages[1])*c1;
        rSecantTensor(0,1) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c2;
        rSecantTensor(1,0) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c2;
        rSecantTensor(2,2) = std::sqrt((1.0 - rDamages[0])*(1.0 - rDamages[1]))*c3;
    }
}


/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<6>>>>;

template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<3>>>>;

} // namespace Kratos
