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
#include "custom_constitutive/generic_small_strain_d_plus_d_minus_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

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

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{

    // Integrate Stress Damage
    Vector& integrated_stress_vector = rValues.GetStressVector();
    array_1d<double, VoigtSize> auxiliar_integrated_stress_vector = integrated_stress_vector;
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
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
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        // S0 = C:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    
        double& tension_threshold = this->GetTensionThreshold();
        double& tension_damage    = this->GetTensionDamage();
        double& compression_threshold = this->GetCompressionThreshold();
        double& compression_damage    = this->GetCompressionDamage();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    //ProcessInfo dummy_process_info;
    //ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    //// We call the integrator
    //double initial_threshold;
    //TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    //this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
	/*this->SetDamage(this->GetNonConvDamage());
	this->SetThreshold(this->GetNonConvThreshold());*/
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
	FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
bool GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == THRESHOLD) {
        return true;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
bool GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    Has(const Variable<Vector>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
bool GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // todo
    //if (rThisVariable == DAMAGE) {
    //    mDamage = rValue;
    //} else if (rThisVariable == THRESHOLD) {
    //    mThreshold = rValue;
    //} else if (rThisVariable == UNIAXIAL_STRESS) {
    //    mUniaxialStress = rValue;
    //} else {
    //    return BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    //}
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
double& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    //if (rThisVariable == DAMAGE) {
    //    rValue = mDamage;
    //} else if (rThisVariable == THRESHOLD) {
    //    rValue = mThreshold;
    //} else if (rThisVariable == UNIAXIAL_STRESS) {
    //    rValue = mUniaxialStress;
    //} else {
    //    return BaseType::GetValue(rThisVariable, rValue);
    //}

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
Vector& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
Matrix& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
double& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
Vector& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
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

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
Matrix& GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    //if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
    //    //1.-Compute total deformation gradient
    //    const Matrix& deformation_gradient_F = rParameterValues.GetDeformationGradientF();
    //    //2.-Right Cauchy-Green tensor C
    //    Matrix right_cauchy_green = prod(trans(deformation_gradient_F), deformation_gradient_F);
    //    Vector strain_vector = ZeroVector(6);

    //    //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
    //    strain_vector[0] = 0.5 * (right_cauchy_green(0, 0) - 1.00);
    //    strain_vector[1] = 0.5 * (right_cauchy_green(1, 1) - 1.00);
    //    strain_vector[2] = 0.5 * (right_cauchy_green(2, 2) - 1.00);
    //    strain_vector[3] = right_cauchy_green(0, 1); // xy
    //    strain_vector[4] = right_cauchy_green(1, 2); // yz
    //    strain_vector[5] = right_cauchy_green(0, 2); // xz

    //    Matrix constitutive_matrix;
    //    this->CalculateElasticMatrix(constitutive_matrix, rParameterValues);

    //    Vector stress = prod(constitutive_matrix, strain_vector);
    //    stress *= (1.0 - mDamage);
    //    rValue =  MathUtils<double>::StressVectorToTensor(stress);
    //    return rValue;
    //} else if (this->Has(rThisVariable)) {
    //    return this->GetValue(rThisVariable, rValue);
    //} else {
    //    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    //}
    
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
int GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const int check_integrator_tension = TConstLawIntegratorTensionType::Check(rMaterialProperties);
    const int check_integrator_compression = TConstLawIntegratorCompressionType::Check(rMaterialProperties);


    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    
    if ((check_base + check_integrator_tension + check_integrator_compression) > 0) return 1;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/


template class GenericSmallStrainDplusDminusDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>,
    GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;






} // Namespace Kratos







