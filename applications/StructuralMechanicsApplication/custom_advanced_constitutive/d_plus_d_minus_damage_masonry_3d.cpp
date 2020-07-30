// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//                   Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

#include "includes/checks.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_advanced_constitutive/d_plus_d_minus_damage_masonry_3d.h"

namespace Kratos
{

DamageDPlusDMinusMasonry3DLaw::DamageDPlusDMinusMasonry3DLaw() {

}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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
        noalias(damage_parameters.TensionStressVector)     = predictive_stress_vector_tension;
        noalias(damage_parameters.CompressionStressVector) = predictive_stress_vector_compression;

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
bool DamageDPlusDMinusMasonry3DLaw::IntegrateStressTensionIfNecessary(
    const double F_tension,
    DamageParameters& rParameters,
    array_1d<double,VoigtSize>& rIntegratedStressVectorTension,
    array_1d<double,VoigtSize> effective_stress_vector,
    ConstitutiveLaw::Parameters& rValues)
{
    bool is_damaging = false;
    const Flags& r_constitutive_law_options = rValues.GetOptions();
    if (F_tension <= tolerance) { //Elastic Case
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
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
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
bool DamageDPlusDMinusMasonry3DLaw::IntegrateStressCompressionIfNecessary(
    const double F_compression,
    DamageParameters& rParameters,
    array_1d<double,VoigtSize>& rIntegratedStressVectorCompression,
    array_1d<double,VoigtSize> effective_stress_vector,
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
void DamageDPlusDMinusMasonry3DLaw::CalculateIntegratedStressVector(
    Vector& rIntegratedStressVector,
    const DamageParameters& rParameters,
    ConstitutiveLaw::Parameters& rValues)
{
    rIntegratedStressVector = (1.0 - rParameters.DamageTension) * rParameters.TensionStressVector +
                              (1.0 - rParameters.DamageCompression) * rParameters.CompressionStressVector;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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
void DamageDPlusDMinusMasonry3DLaw::CalculateSecantTensor(ConstitutiveLaw::Parameters& rValues, Matrix& rSecantTensor)
{
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, rSecantTensor);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::InitializeMaterial(
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
void DamageDPlusDMinusMasonry3DLaw::FinalizeSolutionStep(
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
void DamageDPlusDMinusMasonry3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry3DLaw::Has(const Variable<double>& rThisVariable)
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
bool DamageDPlusDMinusMasonry3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        return true;
    }
    return BaseType::Has(rThisVariable);
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry3DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    if(rThisVariable == INTEGRATED_STRESS_TENSOR) {
        // explicitly returning "false", so we know we must call CalculateValue(...)
        return false;
    }
    return BaseType::Has(rThisVariable);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::SetValue(
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
void DamageDPlusDMinusMasonry3DLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mTensionDamage = rValue[0];
        mTensionThreshold = rValue[1];
        mCompressionDamage = rValue[2];
        mCompressionThreshold = rValue[3];
        mCompressionUniaxialStress = rValue[4];
        mTensionUniaxialStress = rValue[5];
    }
}
/***********************************************************************************/
/***********************************************************************************/
double& DamageDPlusDMinusMasonry3DLaw::GetValue(
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
Vector& DamageDPlusDMinusMasonry3DLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue.resize(6);
        rValue[0] = mTensionDamage;
        rValue[1] = mTensionThreshold;
        rValue[2] = mCompressionDamage;
        rValue[3] = mCompressionThreshold;
        rValue[4] = mCompressionUniaxialStress;
        rValue[5] = mTensionUniaxialStress;
    }
    return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
Matrix& DamageDPlusDMinusMasonry3DLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return BaseType::GetValue(rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
double& DamageDPlusDMinusMasonry3DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return this->GetValue(rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
Vector& DamageDPlusDMinusMasonry3DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/
Matrix& DamageDPlusDMinusMasonry3DLaw::CalculateValue(
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
int DamageDPlusDMinusMasonry3DLaw::Check(
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
    KRATOS_CHECK_VARIABLE_KEY(TRIAXIAL_COMPRESSION_COEFFICIENT);

    if (check_base > 0) return 1;
    return 0;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateEquivalentStressTension(
    array_1d<double, VoigtSize>& rPredictiveStressVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
    const double yield_compression = r_material_properties[YIELD_STRESS_COMPRESSION];
    const double biaxial_compression_multiplier = r_material_properties[BIAXIAL_COMPRESSION_MULTIPLIER];
    const double alpha = (biaxial_compression_multiplier - 1.0)/(2 * biaxial_compression_multiplier - 1.0);
    const double alpha_factor = 1.0 / (1.0 - alpha);
    const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);

    double I1,J2;
    ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
    array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);

    array_1d<double, Dimension> principal_stress_vector;
    ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stress_vector, rPredictiveStressVector);
    const double principal_stress_1 = principal_stress_vector[0];

    if (principal_stress_1 > 0.0){
        rEquivalentStress = alpha_factor * (alpha*I1 + std::sqrt(3.0 * J2) + beta * principal_stress_1) * (yield_tension / yield_compression);
    }
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateEquivalentStressCompression(
    array_1d<double, VoigtSize>& rPredictiveStressVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
    const double yield_compression = r_material_properties[YIELD_STRESS_COMPRESSION];
    const double biaxial_compression_multiplier = r_material_properties[BIAXIAL_COMPRESSION_MULTIPLIER];
    const double shear_compression_reductor = r_material_properties[SHEAR_COMPRESSION_REDUCTOR];
    const double rho = r_material_properties[TRIAXIAL_COMPRESSION_COEFFICIENT];

    KRATOS_ERROR_IF(shear_compression_reductor < 0.0)<< "The SHEAR_COMPRESSION_REDUCTOR is supposed to be a value between 0.0 and 1.0" << std::endl;
    KRATOS_ERROR_IF(shear_compression_reductor > 1.0)<< "The SHEAR_COMPRESSION_REDUCTOR is supposed to be a value between 0.0 and 1.0" << std::endl;
    KRATOS_ERROR_IF(rho <= 0.5)<< "The TRIAXIAL_COMPRESSION_COEFFICIENT is supposed to be a value between 0.5 and 1.0" << std::endl;
    KRATOS_ERROR_IF(rho > 1.0)<< "The TRIAXIAL_COMPRESSION_COEFFICIENT is supposed to be a value between 0.5 and 1.0" << std::endl;

    const double alpha = (biaxial_compression_multiplier - 1.0)/(2.0* biaxial_compression_multiplier - 1.0);
    const double alpha_factor = 1.0 / (1.0 - alpha);
    const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
    const double gamma = 3.0 * (1.0 - rho) / (2.0 * rho - 1.0);

    double I1,J2;
    ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
    array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);

    array_1d<double, Dimension> principal_stress_vector;
    ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stress_vector, rPredictiveStressVector);
    const double principal_stress_1 = principal_stress_vector[0];
    const double principal_stress_3 = principal_stress_vector[2];
    const double smax_macaulay = std::max(principal_stress_1, 0.0);
    const double smax_macaulay_neg = std::abs(std::min(principal_stress_1, 0.0));

    if (principal_stress_3 < 0.0){
        rEquivalentStress = alpha_factor * (alpha*I1 + std::sqrt(3.0 * J2) +
                                            beta * shear_compression_reductor * smax_macaulay +
                                            gamma * smax_macaulay_neg);
    }
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::IntegrateStressVectorTension(
    array_1d<double,VoigtSize>& rPredictiveStressVector,
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
void DamageDPlusDMinusMasonry3DLaw::CalculateDamageParameterTension(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const double Gf = r_material_properties[FRACTURE_ENERGY_TENSION];
    const double E = r_material_properties[YOUNG_MODULUS];
    const double yield_tension = r_material_properties[YIELD_STRESS_TENSION];
    const double l_mat = 2.0 * E * Gf / (std::pow(yield_tension, 2));

    KRATOS_ERROR_IF(CharacteristicLength >= l_mat) << "FRACTURE_ENERGY_TENSION is too low:  2*E*Gt/(ft*ft) = " << l_mat
        << ",   Characteristic Length = " << CharacteristicLength << std::endl;

    rAParameter = 2.0 * (CharacteristicLength / (l_mat - CharacteristicLength));
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::CalculateExponentialDamageTension(
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
void DamageDPlusDMinusMasonry3DLaw::IntegrateStressVectorCompression(
    array_1d<double,VoigtSize>& rPredictiveStressVector,
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
void DamageDPlusDMinusMasonry3DLaw::CalculateBezier3DamageCompression(
    const double UniaxialStress,
    double& rDamage,
    double& rThreshold,
    const double CharacteristicLength,
    ConstitutiveLaw::Parameters& rValues)
{
    // Call the Material Properties
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double young_modulus = r_material_properties[YOUNG_MODULUS];
    const double stress_damage_onset = r_material_properties[DAMAGE_ONSET_STRESS_COMPRESSION];
    const double yield_stress_compression = r_material_properties[YIELD_STRESS_COMPRESSION];
    const double yield_strain_compression = r_material_properties[YIELD_STRAIN_COMPRESSION];
    const double residual_stress_compression = r_material_properties[RESIDUAL_STRESS_COMPRESSION];
    const double bezier_controller_c1 = r_material_properties[BEZIER_CONTROLLER_C1];
    const double bezier_controller_c2 = r_material_properties[BEZIER_CONTROLLER_C2];
    const double bezier_controller_c3 = r_material_properties[BEZIER_CONTROLLER_C3];
    const double fracture_energy_compression = r_material_properties[FRACTURE_ENERGY_COMPRESSION];

    // Calculate missing Bezier Determinators
    const double bezier_control_alpha = 2.0 * (yield_strain_compression - (yield_stress_compression / young_modulus));
    const double strain_damage_onset = stress_damage_onset / young_modulus;
    const double bezier_control_strain_i = yield_stress_compression / young_modulus;
    const double bezier_control_stress_k = residual_stress_compression +
                 (yield_stress_compression - residual_stress_compression) * bezier_controller_c1;
    double bezier_control_strain_j = yield_strain_compression + bezier_control_alpha * bezier_controller_c2;
    double bezier_control_strain_k = 3.0 * yield_strain_compression - 2.0 * yield_stress_compression / young_modulus;
    double bezier_control_strain_r = ( (bezier_control_strain_k - bezier_control_strain_j) *
           (yield_stress_compression - residual_stress_compression)/(yield_stress_compression - bezier_control_stress_k) )
           + bezier_control_strain_j;
    double bezier_control_strain_u = bezier_control_strain_r * bezier_controller_c3;
    const double specific_fracture_energy_compression = fracture_energy_compression / CharacteristicLength;

    // Perform the Energy Regularization of the Bezier Determinators
    this->RegulateBezierDeterminators(
        specific_fracture_energy_compression,
        yield_stress_compression, bezier_control_stress_k, residual_stress_compression, yield_strain_compression,
        bezier_control_strain_j, bezier_control_strain_k, bezier_control_strain_r, bezier_control_strain_u);

    // Compute rDamage
    const double strain_like_counterpart = UniaxialStress / young_modulus;
    double damage_variable_bezier = UniaxialStress;
    if (strain_like_counterpart <= yield_strain_compression) {
        damage_variable_bezier = this->EvaluateBezierCurve(
            strain_like_counterpart,
            strain_damage_onset, bezier_control_strain_i, yield_strain_compression,
            stress_damage_onset, yield_stress_compression, yield_stress_compression);
    } else if (strain_like_counterpart <= bezier_control_strain_k) {
        damage_variable_bezier = this->EvaluateBezierCurve(
            strain_like_counterpart,
            yield_strain_compression, bezier_control_strain_j, bezier_control_strain_k,
            yield_stress_compression, yield_stress_compression, bezier_control_stress_k);
    } else if (strain_like_counterpart <= bezier_control_strain_u) {
        damage_variable_bezier = this->EvaluateBezierCurve(
            strain_like_counterpart,
            bezier_control_strain_k, bezier_control_strain_r, bezier_control_strain_u,
            bezier_control_stress_k, residual_stress_compression, residual_stress_compression);
    } else {
        damage_variable_bezier = residual_stress_compression;
    }

    rDamage = 1.0 - damage_variable_bezier / UniaxialStress;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::RegulateBezierDeterminators(
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

    const double bezier_stretcher = ((specific_dissipated_fracture_energy - bezier_energy_1) /
                                   (BezierEnergy - bezier_energy_1)) - 1.0;

    KRATOS_ERROR_IF(bezier_stretcher <= -1.0) << "Error in Compression Damage: FRACTURE_ENERGY_COMPRESSION is too low, increase it to avoid constitutive snap-back!" << std::endl;

    // Update Strain values
    ej += bezier_stretcher * (ej - ep);
    ek += bezier_stretcher * (ek - ep);
    er += bezier_stretcher * (er - ep);
    eu += bezier_stretcher * (eu - ep);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry3DLaw::ComputeBezierEnergy(
    double& rBezier_energy,
    const double x1, const double x2, const double x3,
    const double y1, const double y2, const double y3)
{
    rBezier_energy = (x2*y1/3.0) + (x3*y1/6.0) - (x2*y3/3) + (x3*y2/3) + (x3*y3/2.0) - x1*((y1/2.0) + (y2/3.0) + (y3/6.0));
}
/***********************************************************************************/
/***********************************************************************************/
double DamageDPlusDMinusMasonry3DLaw::EvaluateBezierCurve(
    const double Xi,
    const double x1, double x2, const double x3,
    const double y1, const double y2, const double y3)
{
    double A = x1 - 2.0 * x2 + x3;
    double B = 2.0 * (x2 - x1);
    double C = x1 - Xi;
    if (std::abs(A) < 1.0e-12) {
        x2 = x2 + 1.0E-6 * (x3-x1);
        A =  x1 - 2.0 * x2 + x3;
        B = 2.0 * (x2 - x1);
        C = x1 - Xi;
    }

    const double D = B * B - 4.0 * A * C;
    const double t = (-B + std::sqrt(D)) / (2.0 * A);
    const double bezier_damage_parameter =  (y1 - 2.0 * y2 + y3) * t * t + (y2 - y1) * 2.0 * t + y1;
    return bezier_damage_parameter;
}
/***********************************************************************************/
/***********************************************************************************/

}// namespace Kratos
