// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner 
//  
//

// System includes

// External includes

// Project includes

#include "includes/checks.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/plane_stress_d_plus_d_minus_damage_masonry_2d.h"

namespace Kratos
{

DamageDPlusDMinusMasonry2DLaw::DamageDPlusDMinusMasonry2DLaw() {
    
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
	this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
	this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
	this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
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
		
		// Compute the equivalent uniaxial Stress in tension and compression
		this->CalculateEquivalentStressTension(predictive_stress_vector, damage_parameters.UniaxialTensionStress, rValues);	
		this->CalculateEquivalentStressCompression(predictive_stress_vector, damage_parameters.UniaxialCompressionStress, rValues);
		
		const double F_tension = damage_parameters.UniaxialTensionStress - damage_parameters.ThresholdTension;
        const double F_compression = damage_parameters.UniaxialCompressionStress - damage_parameters.ThresholdCompression;
        const bool is_damaging_tension = this->IntegrateStressTensionIfNecessary(F_tension, damage_parameters, predictive_stress_vector_tension, predictive_stress_vector, rValues);
		const bool is_damaging_compression = this->IntegrateStressCompressionIfNecessary(F_compression, damage_parameters, predictive_stress_vector_compression, predictive_stress_vector, rValues);
		if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            if (is_damaging_tension || is_damaging_compression) { // Perturbations
                this->CalculateTangentTensor(rValues);
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
bool DamageDPlusDMinusMasonry2DLaw::IntegrateStressTensionIfNecessary(
	const double F_tension,
	DamageParameters& rParameters,
	array_1d<double,3>& rIntegratedStressVectorTension,
    array_1d<double,3> effective_stress_vector,
	ConstitutiveLaw::Parameters& rValues)
{
	bool is_damaging = false;
	const Flags& r_constitutive_law_options = rValues.GetOptions();
	if (F_tension<=tolerance){//Elastic Case
		if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            this->SetNonConvTensionDamage(rParameters.DamageTension);
            this->SetNonConvTensionThreshold(rParameters.ThresholdTension);
        }
		rIntegratedStressVectorTension *= (1.0 - rParameters.DamageTension);
	} else { // Increasing damage...
        const double characteristic_length = ConstitutiveLawUtilities<3>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
        // This routine updates the IntegratedStressVectorTension to verify the yield surf
        this->IntegrateStressVectorTension(
            rIntegratedStressVectorTension, 
            rParameters.UniaxialTensionStress, 
            rParameters.DamageTension, 
            rParameters.ThresholdTension, 
            rValues, characteristic_length);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            this->SetNonConvTensionDamage(rParameters.DamageTension);
            this->SetNonConvTensionThreshold(rParameters.UniaxialTensionStress);
        }
        is_damaging = true;
    }
    // Just for Plotting

    double uniaxial_stress_tension = 0.0;
    this->CalculateEquivalentStressTension(effective_stress_vector, uniaxial_stress_tension, rValues);
    this->SetTensionStress(uniaxial_stress_tension);

    return is_damaging;	
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::IntegrateStressCompressionIfNecessary(
	const double F_compression,
	DamageParameters& rParameters,
	array_1d<double,3>& rIntegratedStressVectorCompression,
    array_1d<double,3> effective_stress_vector,
	ConstitutiveLaw::Parameters& rValues)
{
	bool is_damaging = false;
	const Flags& r_constitutive_law_options = rValues.GetOptions();
    if (F_compression <= tolerance) { // Elastic case
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->SetNonConvCompressionDamage(rParameters.DamageCompression);
            this->SetNonConvCompressionThreshold(rParameters.ThresholdCompression);
        }
        rIntegratedStressVectorCompression *= (1.0 - rParameters.DamageCompression);
    } else { // Increasing damage...
	    const double characteristic_length = ConstitutiveLawUtilities<3>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        // This routine updates the IntegratedStressVectorCompression to verify the yield surf
        this->IntegrateStressVectorCompression(
            rIntegratedStressVectorCompression, 
            rParameters.UniaxialCompressionStress, 
            rParameters.DamageCompression, 
            rParameters.ThresholdCompression, 
            rValues, characteristic_length);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            this->SetNonConvCompressionDamage(rParameters.DamageCompression);
            this->SetNonConvCompressionThreshold(rParameters.UniaxialCompressionStress);
        }
        is_damaging =  true;
	}
    // Just for Plotting
    double uniaxial_stress_compression = 0.0;
    this->CalculateEquivalentStressCompression(effective_stress_vector, uniaxial_stress_compression, rValues);
    this->SetCompressionStress(uniaxial_stress_compression);
    return is_damaging;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateIntegratedStressVector(
	Vector& rIntegratedStressVector,
    const DamageParameters& rParameters,
    ConstitutiveLaw::Parameters& rValues)
{
    rIntegratedStressVector = (1.0 - rParameters.DamageTension) * rParameters.TensionStressVector +
                              (1.0 - rParameters.DamageCompression) * rParameters.CompressionStressVector;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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
void DamageDPlusDMinusMasonry2DLaw::CalculateSecantTensor(ConstitutiveLaw::Parameters& rValues, Matrix& rSecantTensor)
{
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, rSecantTensor);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
{
    // We construct the CL parameters
    double initial_threshold_tension = rMaterialProperties[YIELD_STRESS_TENSION];
	double initial_threshold_compression = rMaterialProperties[DAMAGE_ONSET_STRESS_COMPRESSION];
	
    this->SetTensionThreshold(initial_threshold_tension);
    this->SetCompressionThreshold(initial_threshold_compression);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeSolutionStep(
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
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(const Variable<double>& rThisVariable)
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
bool DamageDPlusDMinusMasonry2DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
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
double& DamageDPlusDMinusMasonry2DLaw::GetValue(
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
Vector& DamageDPlusDMinusMasonry2DLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
Matrix& DamageDPlusDMinusMasonry2DLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
double& DamageDPlusDMinusMasonry2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return this->GetValue(rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
Vector& DamageDPlusDMinusMasonry2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
Matrix& DamageDPlusDMinusMasonry2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        //1.-Compute total deformation gradient
        const Matrix& deformation_gradient_F = rParameterValues.GetDeformationGradientF();
        //2.-Right Cauchy-Green tensor C
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
int DamageDPlusDMinusMasonry2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
	
	KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
	KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
	KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
	KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY_TENSION);
	KRATOS_CHECK_VARIABLE_KEY(DAMAGE_ONSET_STRESS_COMPRESSION);
	KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);
	KRATOS_CHECK_VARIABLE_KEY(YIELD_STRAIN_COMPRESSION);
	KRATOS_CHECK_VARIABLE_KEY(RESIDUAL_STRESS_COMPRESSION);
	KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY_COMPRESSION);
	KRATOS_CHECK_VARIABLE_KEY(BIAXIAL_COMPRESSION_MULTIPLIER);
	KRATOS_CHECK_VARIABLE_KEY(BEZIER_CONTROLLER_C1);
	KRATOS_CHECK_VARIABLE_KEY(BEZIER_CONTROLLER_C2);
	KRATOS_CHECK_VARIABLE_KEY(BEZIER_CONTROLLER_C3);

    if (check_base > 0) return 1;
    return 0;
}
/***********************************************************************************/
/***********************************************************************************/



//////////////////////////////////////////////////////////////
// From here functions could be moved to an utilities.h
//////////////////////////////////////////////////////////////

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateEquivalentStressTension(
	array_1d<double, 3>& rPredictiveStressVector, 
	double& rEquivalentStress,
	ConstitutiveLaw::Parameters& rValues
	)
{
	const Properties& r_material_properties = rValues.GetMaterialProperties();
	
	double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
	double yield_compression = r_material_properties[YIELD_STRESS_COMPRESSION];
	double biaxial_compression_multiplier = r_material_properties[BIAXIAL_COMPRESSION_MULTIPLIER];
	const double alpha = (biaxial_compression_multiplier - 1.0)/(2 * biaxial_compression_multiplier - 1.0);
	const double alpha_factor = 1.0 / (1.0 - alpha);
	const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
	
	
	double I1,J2;
	ConstitutiveLawUtilities<3>::CalculateI1Invariant(rPredictiveStressVector, I1);
	array_1d<double, 3> deviator = ZeroVector(3);
	ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
	
	array_1d<double, 2> rPrincipalStressVector;
	ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(rPrincipalStressVector, rPredictiveStressVector);
	const double principal_stress_1 = rPrincipalStressVector[0];
	const double principal_stress_2 = rPrincipalStressVector[1];
	
	if (principal_stress_1 > 0.0){
		rEquivalentStress = alpha_factor * (alpha*I1 + std::sqrt(3 * J2) + beta * principal_stress_1) * (yield_tension / yield_compression);
	}	
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateEquivalentStressCompression(
	array_1d<double, 3>& rPredictiveStressVector, 
	double& rEquivalentStress,
	ConstitutiveLaw::Parameters& rValues
	)
{
	const Properties& r_material_properties = rValues.GetMaterialProperties();
	
	double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
	double yield_compression = r_material_properties[YIELD_STRESS_COMPRESSION];
	double biaxial_compression_multiplier = r_material_properties[BIAXIAL_COMPRESSION_MULTIPLIER];
	double shear_compression_reductor = r_material_properties[SHEAR_COMPRESSION_REDUCTOR];
	
    //KRATOS_ERROR_IF(shear_compression_reductor < 0.0)<< "The SHEAR_COMPRESSION_REDUCTOR is supposed to be a value between 0.0 and 1.0" << std::endl;
	//KRATOS_ERROR_IF(shear_compression_reductor > 1.0)<< "The SHEAR_COMPRESSION_REDUCTOR is supposed to be a value between 0.0 and 1.0" << std::endl;
	
	const double alpha = (biaxial_compression_multiplier - 1.0)/(2.0* biaxial_compression_multiplier - 1.0);
	const double alpha_factor = 1.0 / (1.0 - alpha);
	const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
	
	double I1,I2,J2;
	ConstitutiveLawUtilities<3>::CalculateI1Invariant(rPredictiveStressVector, I1);
	array_1d<double, 3> deviator = ZeroVector(VoigtSize);
    ConstitutiveLawUtilities<3>::CalculateI2Invariant(rPredictiveStressVector, I2);
    //double J2_temp = (1.0 / 3.0) * std::pow(I1,2) - I2;
	ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
	
	array_1d<double, 2> rPrincipalStressVector;
	ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(rPrincipalStressVector, rPredictiveStressVector);
	const double principal_stress_1 = rPrincipalStressVector[0];
	const double principal_stress_2 = rPrincipalStressVector[1];
	const double smax_macaulay = std::max(principal_stress_1, 0.0);
	
	if (principal_stress_2 < 0.0){
		rEquivalentStress = alpha_factor * (alpha*I1 + std::sqrt(3 * J2) + beta * shear_compression_reductor * smax_macaulay);
	}
	
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::IntegrateStressVectorTension(
	array_1d<double,3>& rPredictiveStressVector,
    const double UniaxialStress,
    double& rDamage,
    double& rThreshold,
    ConstitutiveLaw::Parameters& rValues,
    const double CharacteristicLength)
{
	double damage_parameter;
	this->CalculateDamageParameterTension(rValues, damage_parameter, CharacteristicLength);
	this->CalculateExponentialDamageTension(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
	
	rPredictiveStressVector *= (1.0 - rDamage);
}
/***********************************************************************************/
/***********************************************************************************/		
void DamageDPlusDMinusMasonry2DLaw::CalculateDamageParameterTension(
	ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength)
{
	const Properties& r_material_properties = rValues.GetMaterialProperties();
	
	const double Gf = r_material_properties[FRACTURE_ENERGY_TENSION];
	const double E = r_material_properties[YOUNG_MODULUS];
	const double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
	const double l_mat = 2.0 * E * Gf / (std::pow(yield_tension, 2));
	rAParameter = 2.0 * (CharacteristicLength / (l_mat - CharacteristicLength));
    if(CharacteristicLength >= l_mat)
			{
				std::stringstream ss;
				ss << "FRACTURE_ENERGY_TENSION is too low:  2*E*Gt/(ft*ft) = " << l_mat
					<< ",   Characteristic Length = " << CharacteristicLength << std::endl;
				std::cout << ss.str();
				exit(-1);
			}
	//KRATOS_ERROR_IF(rAParameter < 0.0) << "FRACTURE_ENERGY_TENSION is too low, increase it ..." << std::endl;	
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateExponentialDamageTension(
	const double UniaxialStress,
	const double Threshold,
	const double DamageParameter,
	const double CharacteristicLength,
	ConstitutiveLaw::Parameters& rValues,
	double& rDamage)
{
	const Properties& r_material_properties = rValues.GetMaterialProperties();
	const double initial_threshold = r_material_properties[YIELD_STRESS_TENSION];
	rDamage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter * (1.0 - (UniaxialStress / initial_threshold)));
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::IntegrateStressVectorCompression(
	array_1d<double,3>& rPredictiveStressVector,
	const double UniaxialStress,
    double& rDamage,
    double& rThreshold,
    ConstitutiveLaw::Parameters& rValues,
    const double CharacteristicLength)
{
	this->CalculateBezier3DamageCompression(UniaxialStress, rDamage, rThreshold, CharacteristicLength, rValues);
	rPredictiveStressVector *= (1.0 - rDamage);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateBezier3DamageCompression(
	const double UniaxialStress,
	double& rDamage,
	double& rThreshold,
	const double CharacteristicLength,
	ConstitutiveLaw::Parameters& rValues)
{
	// Call the Material Properties 
	const Properties& r_material_properties = rValues.GetMaterialProperties();
	const double E = r_material_properties[YOUNG_MODULUS]; 
	const double s0 = r_material_properties[DAMAGE_ONSET_STRESS_COMPRESSION];
	const double sp = r_material_properties[YIELD_STRESS_COMPRESSION];
	const double ep = r_material_properties[YIELD_STRAIN_COMPRESSION];
	const double sr = r_material_properties[RESIDUAL_STRESS_COMPRESSION];
	const double c1 = r_material_properties[BEZIER_CONTROLLER_C1];
	const double c2 = r_material_properties[BEZIER_CONTROLLER_C2];
	const double c3 = r_material_properties[BEZIER_CONTROLLER_C3];
	const double Gc = r_material_properties[FRACTURE_ENERGY_COMPRESSION];
	
	// Calculate missing Bezier Determinators
	const double alpha = 2.0 * (ep - (sp/E));
	const double e0 = s0 / E;
	const double ei = sp / E;
	const double sk = sr + (sp - sr) * c1;
	double ej = ep + alpha * c2;
	//double ek = ej + alpha * (1.0 - c2);
	double ek = 3.0 * ep - 2.0 * sp / E;
	double er = ( (ek - ej) * (sp - sr)/(sp -sk) ) + ej;
	double eu = er * c3;
	const double gc = Gc / CharacteristicLength;
	
	// Perform the Energy Regularization of the Bezier Determinators
	//double gc_bezier; 
	this->RegulateBezierDeterminators(gc, sp, sk, sr, ep, ej, ek, er, eu);
	
	// Compute rDamage
	double StrainLikeCounterpart = UniaxialStress / E;
	double DamageVariableBezier = UniaxialStress;
	if (StrainLikeCounterpart <= ep){
		DamageVariableBezier = this->EvaluateBezierCurve(StrainLikeCounterpart, e0, ei, ep, s0, sp, sp);
	} else {if (StrainLikeCounterpart <= ek){
			DamageVariableBezier = this->EvaluateBezierCurve(StrainLikeCounterpart, ep, ej, ek, sp, sp, sk);
			} else {if (StrainLikeCounterpart <= eu){
				DamageVariableBezier = this->EvaluateBezierCurve(StrainLikeCounterpart, ek, er, eu, sk, sr, sr);
				} else {
					DamageVariableBezier = sr;
					   }
				   }
		   }
	rDamage = 1.0 - DamageVariableBezier / UniaxialStress;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::RegulateBezierDeterminators(
	const double specific_dissipated_fracture_energy,
	const double sp, const double sk, const double sr, const double ep,
	double& ej, double& ek, double& er, double& eu)
{	
	const double bezier_energy_1 = sp * ep / 2.0;
	double bezier_energy_2;
	double bezier_energy_3;
	this->ComputeBezierEnergy(bezier_energy_2, ep, ej, ek, sp, sp, sk);
	this->ComputeBezierEnergy(bezier_energy_3, ek, er, eu, sk, sr, sr);
	const double BezierEnergy = bezier_energy_1 + bezier_energy_2 + bezier_energy_3;
	
	const double BezierStretcher = ((specific_dissipated_fracture_energy - bezier_energy_1) / 
								   (BezierEnergy - bezier_energy_1)) - 1.0;

	if (BezierStretcher <= -1.0)
	{
		std::stringstream ss;
		ss << "Error in Compression Damage: FRACTURE_ENERGY_COMPRESSION is too low, increase it to avoid constitutive snap-back!" << std::endl;
		ss << "Input Gc/lch = " << specific_dissipated_fracture_energy << std::endl;
		std::cout << ss.str();
		exit(-1);
	}
	//KRATOS_ERROR_IF(BezierStretcher <= -1.0) << "FRACTURE_ENERGY_COMPRESSION is too low, increase it to avoid constitutive snap-back!" << std::endl;	
								   
    // Update Strain values
	ej = ej + BezierStretcher * (ej - ep);
	ek = ek + BezierStretcher * (ek - ep);
	er = er + BezierStretcher * (er - ep);
	eu = eu + BezierStretcher * (eu - ep);

	double CheckBezierEnergy2, CheckBezierEnergy3;
    this-> ComputeBezierEnergy(CheckBezierEnergy2, ep, ej, ek, sp, sp, sk);
    this->ComputeBezierEnergy(CheckBezierEnergy3, ek, er, eu, sk, sr, sr);
    double CheckBezierEnergy = CheckBezierEnergy2 + CheckBezierEnergy3 + bezier_energy_1;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::ComputeBezierEnergy(
	double& BezierG,
	const double x1, const double x2, const double x3,
	const double y1, const double y2, const double y3)
{  
	BezierG = (x2*y1/3.0) + (x3*y1/6.0) - (x2*y3/3) + (x3*y2/3) + (x3*y3/2.0) - x1*((y1/2.0) + (y2/3.0) + (y3/6.0));
}
/***********************************************************************************/
/***********************************************************************************/
double DamageDPlusDMinusMasonry2DLaw::EvaluateBezierCurve(
    const double Xi, 
    const double x1, double x2, const double x3, 
    const double y1, const double y2, const double y3)
{
    double A = x1 - 2.0 * x2 + x3;
    double B = 2.0 * (x2 - x1);
    double C = x1 - Xi;
	double AA = x3 * x1 - std::pow(x2, 2);
	double KK = AA / A;
    if (std::abs(A) < 1.0e-12)
    {
        x2 = x2 + 1.0E-6 * (x3-x1);
        A =  x1 - 2.0 * x2 + x2;
        B = 2.0 * (x2 - x1);
        C = x1 - Xi;
    }
    double D = B * B - 4.0 * A * C;
    double t = (-B + std::sqrt(D)) / (2.0 * A);
    double bezier_damage_parameter =  (y1 - 2.0 * y2 + y3) * t * t + (y2 - y1) * 2.0 * t + y1;
    return bezier_damage_parameter;
}

/***********************************************************************************/
/***********************************************************************************/

}// namespace Kratos