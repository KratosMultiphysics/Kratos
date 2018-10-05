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
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_compression_constitutive_law_integrator.h"
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_tension_constitutive_law_integrator.h"

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
    KRATOS_TRY

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

        DamageParameters damage_parameters;
        damage_parameters.ThresholdTension = this->GetTensionThreshold();
        damage_parameters.DamageTension    = this->GetTensionDamage();
        damage_parameters.ThresholdCompression = this->GetCompressionThreshold();
        damage_parameters.DamageCompression    = this->GetCompressionDamage();

        // S0 = C0:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);

        // Perform the separation of the Stress in tension and compression
        array_1d<double, VoigtSize> predictive_stress_vector_tension, predictive_stress_vector_compression;
        ConstitutiveLawUtilities<VoigtSize>::SpectralDecomposition(predictive_stress_vector, predictive_stress_vector_tension, predictive_stress_vector_compression);

        damage_parameters.TensionStressVector     = predictive_stress_vector_tension;
        damage_parameters.CompressionStressVector = predictive_stress_vector_compression;
        
        TConstLawIntegratorTensionType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector_tension, r_strain_vector,
            damage_parameters.UniaxialTensionStress, rValues);
        TConstLawIntegratorCompressionType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector_compression, r_strain_vector,
            damage_parameters.UniaxialCompressionStress, rValues);

        const double F_tension = damage_parameters.UniaxialTensionStress - damage_parameters.ThresholdTension;
        const double F_compression = damage_parameters.UniaxialCompressionStress - damage_parameters.ThresholdCompression;

        const bool is_damaging_tension = this->IntegrateStressTensionIfNecessary(F_tension, damage_parameters, predictive_stress_vector_tension, rValues);
        const bool is_damaging_compression = this->IntegrateStressCompressionIfNecessary(F_compression, damage_parameters, predictive_stress_vector_compression, rValues);
  
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            if (is_damaging_tension || is_damaging_compression) { // Perturbations
                this->CalculateTangentTensor(rValues);
                noalias(r_tangent_tensor) = rValues.GetConstitutiveMatrix();
            } else { // Secant matrix
                this->CalculateSecantTensor(rValues, r_tangent_tensor);
            }
        }
		this->CalculateIntegratedStressVector(integrated_stress_vector, damage_parameters, rValues);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
bool GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    IntegrateStressTensionIfNecessary(
        const double F_tension, 
        DamageParameters& rParameters, 
        array_1d<double, VoigtSize>& IntegratedStressVectorTension,
        ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    bool is_damaging = false;
    if (F_tension <= 0.0) { // Elastic case
        this->SetNonConvTensionDamage(rParameters.DamageTension);
        this->SetNonConvTensionThreshold(rParameters.ThresholdTension);
        IntegratedStressVectorTension *= (1.0 - rParameters.DamageTension);
    } else { // Increasing damage...
        const double characteristic_length = rValues.GetElementGeometry().Length();

        // This routine updates the IntegratedStressVectorTension to verify the yield surf
        TConstLawIntegratorTensionType::IntegrateStressVector(IntegratedStressVectorTension, rParameters.UniaxialTensionStress, 
            rParameters.DamageTension, rParameters.ThresholdTension, rValues, characteristic_length);
        this->SetNonConvTensionDamage(rParameters.DamageTension);
        this->SetNonConvTensionThreshold(rParameters.UniaxialTensionStress);
        is_damaging =  true;
    }

    // Just for plotting...
    double uniaxial_stress_tension = 0.0;
    TConstLawIntegratorTensionType::YieldSurfaceType::CalculateEquivalentStress(IntegratedStressVectorTension, rValues.GetStrainVector(),
                                                                                uniaxial_stress_tension, rValues);
    uniaxial_stress_tension /= TConstLawIntegratorTensionType::YieldSurfaceType::GetScaleFactorTension(rValues.GetMaterialProperties());
    this->SetTensionStress(uniaxial_stress_tension);

    return is_damaging;
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
bool GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    IntegrateStressCompressionIfNecessary(
        const double F_compression, 
        DamageParameters& rParameters, 
        array_1d<double, VoigtSize>& IntegratedStressVectorCompression,
        ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    bool is_damaging = false;
    if (F_compression <= 0.0) { // Elastic case
        this->SetNonConvCompressionDamage(rParameters.DamageCompression);
        this->SetNonConvCompressionThreshold(rParameters.ThresholdCompression);
        IntegratedStressVectorCompression *= (1.0 - rParameters.DamageCompression);
    } else { // Increasing damage...
        const double characteristic_length = rValues.GetElementGeometry().Length();

        // This routine updates the IntegratedStressVectorCompression to verify the yield surf
        TConstLawIntegratorCompressionType::IntegrateStressVector(IntegratedStressVectorCompression, rParameters.UniaxialCompressionStress, 
            rParameters.DamageCompression, rParameters.ThresholdCompression, rValues, characteristic_length);
        this->SetNonConvCompressionDamage(rParameters.DamageCompression);
        this->SetNonConvCompressionThreshold(rParameters.UniaxialCompressionStress);
        is_damaging =  true;
    }

    // Just for plotting...
    double uniaxial_stress_compression = 0.0;
    TConstLawIntegratorCompressionType::YieldSurfaceType::CalculateEquivalentStress(IntegratedStressVectorCompression, rValues.GetStrainVector(),
                                                                                uniaxial_stress_compression, rValues);
    this->SetCompressionStress(uniaxial_stress_compression);

    return is_damaging;
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    CalculateIntegratedStressVector(
        Vector& rIntegratedStressVector,
        const DamageParameters& rParameters,
        ConstitutiveLaw::Parameters& rValues)
{
    // Matrix projection_operator, A;
    // const Vector strain_vector = rValues.GetStrainVector();
    // const Matrix strain_tensor = MathUtils<double>::StrainVectorToTensor(strain_vector);
    // ConstitutiveLawUtilities<VoigtSize>::CalculateProjectionOperator(strain_vector, projection_operator);
    // Matrix identity_matrix = IdentityMatrix(Dimension, Dimension);

    // A = std::sqrt(1.0 - rParameters.DamageTension) * projection_operator + std::sqrt(1.0 - rParameters.DamageCompression) * (identity_matrix - projection_operator);

    // const Properties material_props = rValues.GetMaterialProperties();
    // const double E  = material_props[YOUNG_MODULUS];
    // const double nu = material_props[POISSON_RATIO];
    // Matrix constitutive_tensor = ZeroMatrix(Dimension, Dimension);
	// this->CalculateLinearElasticMatrix(constitutive_tensor, E, nu);

    // Provisional!!!! TODO-> convertir A en voigt y S = A*D*A*E
    rIntegratedStressVector = (1.0 - rParameters.DamageTension) * rParameters.TensionStressVector +
                              (1.0 - rParameters.DamageCompression) * rParameters.CompressionStressVector;
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
    CalculateSecantTensor(ConstitutiveLaw::Parameters& rValues, Matrix& rSecantTensor)
{
    // TODO improve this
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, rSecantTensor);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
void GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>::
    InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrators
    double initial_threshold_tension, initial_threshold_compression;
    TConstLawIntegratorTensionType::GetInitialUniaxialThreshold(aux_param, initial_threshold_tension);
    this->SetTensionThreshold(initial_threshold_tension);
    
    TConstLawIntegratorCompressionType::GetInitialUniaxialThreshold(aux_param, initial_threshold_compression);
    this->SetCompressionThreshold(initial_threshold_compression);
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
    this->SetTensionDamage(this->GetNonConvTensionDamage());
    this->SetTensionThreshold(this->GetNonConvTensionThreshold());

    this->SetCompressionDamage(this->GetNonConvCompressionDamage());
    this->SetCompressionThreshold(this->GetNonConvCompressionThreshold());
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
    if (rThisVariable == DAMAGE_TENSION) {
        return true;
    } else if (rThisVariable == THRESHOLD_TENSION) {
        return true;
    } else if (rThisVariable == DAMAGE_COMPRESSION) {
		return true;
    } else if (rThisVariable == THRESHOLD_COMPRESSION) {
		return true;
    } else if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION) {
		return true; 
    } else if (rThisVariable == UNIAXIAL_STRESS_TENSION) {
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
    if (rThisVariable == DAMAGE_TENSION) {
       mTensionDamage = rValue;
    } else if (rThisVariable == THRESHOLD_TENSION) {
       mTensionThreshold = rValue;
    } else if (rThisVariable == DAMAGE_COMPRESSION) {
       mCompressionDamage = rValue;
    } else if (rThisVariable == THRESHOLD_COMPRESSION) {
       mCompressionThreshold = rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION) {
       mCompressionUniaxialStress = rValue; 
    } else if (rThisVariable == UNIAXIAL_STRESS_TENSION) {
       mTensionUniaxialStress = rValue; 
    } else {
       return BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
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
    if (rThisVariable == DAMAGE_TENSION) {
       rValue = mTensionDamage;
    } else if (rThisVariable == THRESHOLD_TENSION) {
       rValue = mTensionThreshold;
    } else if (rThisVariable == DAMAGE_COMPRESSION) {
       rValue = mCompressionDamage;
    } else if (rThisVariable == THRESHOLD_COMPRESSION) {
       rValue = mCompressionThreshold;
    } else if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION) {
       rValue = mCompressionUniaxialStress;
    } else if (rThisVariable == UNIAXIAL_STRESS_TENSION) {
       rValue = mTensionUniaxialStress;
    } else {
       return BaseType::GetValue(rThisVariable, rValue);
    }
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
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        //1.-Compute total deformation gradient
        const Matrix& deformation_gradient_F = rParameterValues.GetDeformationGradientF();
        //2.-Right Cauchy-Green tensor C
		KRATOS_WATCH(deformation_gradient_F)
        Matrix right_cauchy_green = prod(trans(deformation_gradient_F), deformation_gradient_F);
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

        Vector stress = prod(constitutive_matrix, strain_vector);
        //stress *= (1.0 - mDamage);
        rValue =  MathUtils<double>::StressVectorToTensor(stress);
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

template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;

} // Namespace Kratos