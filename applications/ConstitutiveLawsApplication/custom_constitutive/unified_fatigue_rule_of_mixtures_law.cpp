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
//  Main authors:    Sergio Jimenez
//                   Alejandro Cornejo
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "unified_fatigue_rule_of_mixtures_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"


namespace Kratos
{
ConstitutiveLaw::Pointer UnifiedFatigueRuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{
    const double high_cycle_fatigue_initial_volumetric_participation = NewParameters["combination_factors"][0].GetDouble();
    const int voigt_size = 6;
    Vector parallel_directions(voigt_size);
    return Kratos::make_shared<UnifiedFatigueRuleOfMixturesLaw>(high_cycle_fatigue_initial_volumetric_participation, parallel_directions);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient(dimension, dimension);
        this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
        const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
        // Doing resize in case is needed
        if (r_strain_vector.size() != voigt_size)
            r_strain_vector.resize(voigt_size);

         // Identity matrix
        Matrix identity_matrix(dimension, dimension);
        for (IndexType i = 0; i < dimension; ++i) {
            for (IndexType j = 0; j < dimension; ++j) {
                if (i == j) identity_matrix(i, j) = 1.0;
                else identity_matrix(i, j) = 0.0;
            }
        }

        // Calculating the inverse of the left Cauchy tensor
        Matrix inverse_B_tensor(dimension, dimension);
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector hich is equal to component CLs strains
        Vector& r_strain_vector = rValues.GetStrainVector();

        // This method integrates the stress according to each simple material CL
        Vector high_cycle_fatigue_stress_vector, ultra_low_cycle_fatigue_stress_vector;
        this->IntegrateStressesOfHCFAndULCFModels(rValues, r_strain_vector, r_strain_vector, high_cycle_fatigue_stress_vector, ultra_low_cycle_fatigue_stress_vector);

        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mHCFVolumetricParticipation * high_cycle_fatigue_stress_vector
                                     + (1.0 - mHCFVolumetricParticipation) * ultra_low_cycle_fatigue_stress_vector;
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/
void UnifiedFatigueRuleOfMixturesLaw::IntegrateStressesOfHCFAndULCFModels(
    ConstitutiveLaw::Parameters& rValues,
    Vector rHCFStrainVector,
    Vector rULCFStrainVector,
    Vector& rHCFStressVector,
    Vector& rULCFStressVector
)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_HCF_cl = *(it_cl_begin);
    const auto& r_props_ULCF_cl = *(it_cl_begin + 1);

    ConstitutiveLaw::Parameters values_HCF  = rValues;
    ConstitutiveLaw::Parameters values_ULCF = rValues;

    values_HCF.SetStrainVector(rHCFStrainVector);
    values_ULCF.SetStrainVector(rULCFStrainVector);

    // Integrate Stress of the HCF part
    values_HCF.SetMaterialProperties(r_props_HCF_cl);
    mpHCFConstitutiveLaw->CalculateMaterialResponseCauchy(values_HCF);
    rHCFStressVector = values_HCF.GetStressVector();

    // Integrate Stress of the UÑCF part
    values_ULCF.SetMaterialProperties(r_props_ULCF_cl);
    mpULCFConstitutiveLaw->CalculateMaterialResponseCauchy(values_ULCF);
    rULCFStressVector = values_ULCF.GetStressVector();
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& r_strain_vector = rValues.GetStrainVector();
    // mPreviousStrainVector = r_strain_vector;

    // Recalculation to obtain the serial_strain_matrix and store the value
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        // We call the FinalizeMaterialResponse of the HCF and ULCF CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_HCF_cl = *(it_cl_begin);
        const auto& r_props_ULCF_cl = *(it_cl_begin + 1);

        ConstitutiveLaw::Parameters values_HCF  = rValues;
        ConstitutiveLaw::Parameters values_ULCF = rValues;

        values_HCF.SetMaterialProperties(r_props_HCF_cl);
        values_ULCF.SetMaterialProperties(r_props_ULCF_cl);

        values_HCF.SetStrainVector(r_strain_vector);
        values_ULCF.SetStrainVector(r_strain_vector);

        mpHCFConstitutiveLaw->FinalizeMaterialResponseCauchy(values_HCF);
        mpULCFConstitutiveLaw->FinalizeMaterialResponseCauchy(values_ULCF);
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& UnifiedFatigueRuleOfMixturesLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return mpHCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return mpULCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Vector& UnifiedFatigueRuleOfMixturesLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return mpHCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return mpULCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& UnifiedFatigueRuleOfMixturesLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return mpHCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return mpULCFConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool UnifiedFatigueRuleOfMixturesLaw::Has(const Variable<bool>& rThisVariable)
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool UnifiedFatigueRuleOfMixturesLaw::Has(const Variable<double>& rThisVariable)
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool UnifiedFatigueRuleOfMixturesLaw::Has(const Variable<Vector>& rThisVariable)
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool UnifiedFatigueRuleOfMixturesLaw::Has(const Variable<Matrix>& rThisVariable)
{
    if (mpHCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpULCFConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& UnifiedFatigueRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if (rThisVariable == UNIAXIAL_STRESS_HCF) {
        return mpHCFConstitutiveLaw->CalculateValue(rParameterValues, UNIAXIAL_STRESS, rValue);
    } else if (rThisVariable == UNIAXIAL_STRESS_ULCF) {
        return mpULCFConstitutiveLaw->CalculateValue(rParameterValues, UNIAXIAL_STRESS, rValue);
    } else {
        return this->GetValue(rThisVariable, rValue);
    }

}

/***********************************************************************************/
/***********************************************************************************/

Vector& UnifiedFatigueRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_HCF_cl = *(it_cl_begin);
    const auto r_props_ULCF_cl  = *(it_cl_begin + 1);

    KRATOS_ERROR_IF_NOT(r_props_HCF_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    KRATOS_ERROR_IF_NOT(r_props_ULCF_cl.Has(CONSTITUTIVE_LAW))  << "No constitutive law set" << std::endl;

    mpHCFConstitutiveLaw = r_props_HCF_cl[CONSTITUTIVE_LAW]->Clone();
    mpULCFConstitutiveLaw  = r_props_ULCF_cl[CONSTITUTIVE_LAW]->Clone();
    mpHCFConstitutiveLaw->InitializeMaterial(r_props_HCF_cl, rElementGeometry, rShapeFunctionsValues);
    mpULCFConstitutiveLaw ->InitializeMaterial(r_props_ULCF_cl, rElementGeometry, rShapeFunctionsValues);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& UnifiedFatigueRuleOfMixturesLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We do some special operations for constitutive matrices
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

         // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == DEFORMATION_GRADIENT) { // TODO: Make in the future modifications for take into account different layers combinations
        noalias(rValue) = rParameterValues.GetDeformationGradientF();
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseCauchy(rParameterValues);
        rValue = MathUtils<double>::StressVectorToTensor(rParameterValues.GetStressVector());

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;
    } else {
        Matrix aux_value;
        Properties material_properties  = rParameterValues.GetMaterialProperties();
        Properties& r_prop = material_properties.GetSubProperties(0);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpHCFConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mHCFVolumetricParticipation) * aux_value;

        r_prop = material_properties.GetSubProperties(1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpHCFConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mHCFVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }
    return(rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueRuleOfMixturesLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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
} // namespace Kratos
