// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:                 BSD License
//                               license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo Velazquez
//                   Sergio Jimenez Reyes
//                   Lucia Gratiela Barbu
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/composites/traction_separation_law.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D()
    : BaseType()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const std::vector<double>& rCombinationFactors) : BaseType(rCombinationFactors)
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const TractionSeparationLaw3D<TDim>& rOther)
    : BaseType(rOther),
      mDelaminationDamageModeOne(rOther.mDelaminationDamageModeOne),
      mDelaminationDamageModeTwo(rOther.mDelaminationDamageModeTwo),
      mThresholdModeOne(rOther.mThresholdModeOne),
      mThresholdModeTwo(rOther.mThresholdModeTwo)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::Pointer TractionSeparationLaw3D<TDim>::Clone() const
{
    return Kratos::make_shared<TractionSeparationLaw3D>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::Pointer TractionSeparationLaw3D<TDim>::Create(Kratos::Parameters NewParameters) const
{
    // We do some checks
    KRATOS_ERROR_IF_NOT(NewParameters.Has("combination_factors")) << "TractionSeparationLaw3D: Please define combination_factors" << std::endl;

    const SizeType number_of_factors = NewParameters["combination_factors"].size();

    // We create the vectors
    std::vector<double> combination_factors(number_of_factors);

    for (IndexType i_layer = 0; i_layer < number_of_factors; ++i_layer) {
        combination_factors[i_layer] = NewParameters["combination_factors"][i_layer].GetDouble();
    }

    KRATOS_ERROR_IF(number_of_factors == 0) << "Please define the combination factors" << std::endl;

    // We create the law
    return Kratos::make_shared<TractionSeparationLaw3D>(combination_factors);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::~TractionSeparationLaw3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    bool has = false;

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_ONE) {
        has = true;
    } else if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_TWO) {
        has = true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    const auto& r_combination_factors = this->GetCombinationFactors();

    rValue.clear();

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_ONE) {

        rValue.resize(r_combination_factors.size()+1, false);

        noalias(rValue) = mDelaminationDamageModeOne;
        return rValue;
    } else if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_TWO) {

        rValue.resize(r_combination_factors.size()+1, false);

        noalias(rValue) = mDelaminationDamageModeTwo;
        return rValue;
    } else {
        BaseType::GetValue(rThisVariable,rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it layer by layer
    bool valid_input = true;
    const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
    for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
        if (p_law->ValidateInput(r_prop)) {
            valid_input = false;
            break;
        }
    }

    return valid_input;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const ConstitutiveLaw::GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();

    BaseType::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

    mDelaminationDamageModeOne.resize(r_p_constitutive_law_vector.size()+1, false);
    noalias(mDelaminationDamageModeOne) = ZeroVector(r_p_constitutive_law_vector.size()+1);

    mDelaminationDamageModeTwo.resize(r_p_constitutive_law_vector.size()+1, false);
    noalias(mDelaminationDamageModeTwo) = ZeroVector(r_p_constitutive_law_vector.size()+1);

    mThresholdModeOne.resize(r_p_constitutive_law_vector.size()-1, false);
    for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        mThresholdModeOne[i] = rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH];
    }

    mThresholdModeTwo.resize(r_p_constitutive_law_vector.size()-1, false);
    for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        mThresholdModeTwo[i] = rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is(BaseType::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();
        Vector auxiliar_stress_vector  = ZeroVector(VoigtSize);
        Vector delamination_damage_affected_stress_vector  = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
        const auto& r_combination_factors = this->GetCombinationFactors();

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }

        std::vector<bool> negative_interfacial_stress_indicator(r_p_constitutive_law_vector.size()+1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()+1; ++i) {
            negative_interfacial_stress_indicator[i] = false;
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        const double tolerance = std::numeric_limits<double>::epsilon();
        Vector DelaminationDamageModeOne(r_p_constitutive_law_vector.size()+1);
        Vector DelaminationDamageModeTwo(r_p_constitutive_law_vector.size()+1);
        Vector ThresholdModeOne(r_p_constitutive_law_vector.size()-1);
        Vector ThresholdModeTwo(r_p_constitutive_law_vector.size()-1);

        noalias(DelaminationDamageModeOne) = mDelaminationDamageModeOne;
        noalias(DelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
        noalias(ThresholdModeOne) = mThresholdModeOne;
        noalias(ThresholdModeTwo) = mThresholdModeTwo;

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {

            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

            const double equivalent_stress_mode_one = interfacial_stress[i][0];
            const double equivalent_stress_mode_two = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));

            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < tolerance) {
                negative_interfacial_stress_indicator[i+1] = true;
            }

            // Damage calculation

            const double T0n = r_material_properties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Ei = r_material_properties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? r_material_properties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? r_material_properties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

            const double F_mode_one = equivalent_stress_mode_one - ThresholdModeOne[i];
            if (F_mode_one > tolerance) {

                DelaminationDamageModeOne[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIc, Ei, T0n, equivalent_stress_mode_one);
            }

            const double F_mode_two = equivalent_stress_mode_two - ThresholdModeTwo[i];
            if (F_mode_two > tolerance) {

                DelaminationDamageModeTwo[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIIc, Gi, T0s, equivalent_stress_mode_two);
            }

            // End damage calculation
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            double layer_damage_variable_mode_one = 0.0;
            double layer_damage_variable_mode_two = 0.0;

            if (DelaminationDamageModeOne[i+1] > DelaminationDamageModeOne[i]) {
                layer_damage_variable_mode_one = DelaminationDamageModeOne[i+1];
            } else {
                layer_damage_variable_mode_one = DelaminationDamageModeOne[i];
            }

            if (DelaminationDamageModeTwo[i+1] > DelaminationDamageModeTwo[i]) {
                layer_damage_variable_mode_two = DelaminationDamageModeTwo[i+1];
            } else {
                layer_damage_variable_mode_two = DelaminationDamageModeTwo[i];
            }

            layer_stress[i][2] *= ((1.0-layer_damage_variable_mode_one));
            layer_stress[i][4] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
            layer_stress[i][5] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            const double factor = r_combination_factors[i];
            delamination_damage_affected_stress_vector += factor * layer_stress[i];
        }

        auxiliar_stress_vector = delamination_damage_affected_stress_vector;

        noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, BaseType::StressMeasure_PK2);
        }

        // Previous flags restored
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is(BaseType::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }

        std::vector<bool> negative_interfacial_stress_indicator(r_p_constitutive_law_vector.size()+1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()+1; ++i) {
            negative_interfacial_stress_indicator[i] = false;
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            p_law->FinalizeMaterialResponsePK2(rValues);

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        const double tolerance = std::numeric_limits<double>::epsilon();
        Vector DelaminationDamageModeOne(r_p_constitutive_law_vector.size()+1);
        Vector DelaminationDamageModeTwo(r_p_constitutive_law_vector.size()+1);
        Vector ThresholdModeOne(r_p_constitutive_law_vector.size()-1);
        Vector ThresholdModeTwo(r_p_constitutive_law_vector.size()-1);

        noalias(DelaminationDamageModeOne) = mDelaminationDamageModeOne;
        noalias(DelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
        noalias(ThresholdModeOne) = mThresholdModeOne;
        noalias(ThresholdModeTwo) = mThresholdModeTwo;

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

            const double equivalent_stress_mode_one = interfacial_stress[i][0];
            const double equivalent_stress_mode_two = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));

            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < tolerance) {
                negative_interfacial_stress_indicator[i+1] = true;
            }

            // Damage calculation
            const double T0n = r_material_properties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Ei = r_material_properties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? r_material_properties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? r_material_properties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

            const double F_mode_one = equivalent_stress_mode_one - ThresholdModeOne[i];
            if (F_mode_one > tolerance) {

                DelaminationDamageModeOne[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIc, Ei, T0n, equivalent_stress_mode_one);

                mDelaminationDamageModeOne[i+1] = DelaminationDamageModeOne[i+1];
                mThresholdModeOne[i] = equivalent_stress_mode_one;
            }

            const double F_mode_two = equivalent_stress_mode_two - ThresholdModeTwo[i];
            if (F_mode_two > tolerance) {

                DelaminationDamageModeTwo[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIIc, Gi, T0s, equivalent_stress_mode_two);

                mDelaminationDamageModeTwo[i+1] = DelaminationDamageModeTwo[i+1];
                mThresholdModeTwo[i] = equivalent_stress_mode_two;
            }

            // End damage calculation
        }

        // Previous flags restored
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double TractionSeparationLaw3D<TDim>::CalculateDelaminationDamageExponentialSoftening(
    ConstitutiveLaw::Parameters& rValues,
    const double GI,
    const double E,
    const double T0,
    const double equivalent_stress)
{
    const double characteristic_length = 0.6343 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry()));
    const double AParameter = 1.00 / (GI * E / (characteristic_length * std::pow(T0, 2)) - 0.5); // Exponential

    KRATOS_ERROR_IF(AParameter < 0.0) << "AParameter is negative." << std::endl;

    double DelaminationDamage = 1.0 - (T0 / equivalent_stress) * std::exp(AParameter * (1.0 - equivalent_stress / T0)); // Exponential

    DelaminationDamage = (DelaminationDamage >= 0.99999) ? 0.99999 : DelaminationDamage;
    DelaminationDamage = (DelaminationDamage < 0.0) ? 0.0 : DelaminationDamage;
    return DelaminationDamage;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int TractionSeparationLaw3D<TDim>::Check(
    const Properties& rMaterialProperties,
    const ConstitutiveLaw::GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    BaseType::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);

    // Check if input parameters are completely defined
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH) || rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR)) << "INTERFACIAL_NORMAL_STRENGTH is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH) || rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR)) << "INTERFACIAL_SHEAR_STRENGTH is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY) || rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) << "MODE_ONE_FRACTURE_ENERGY is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY) || rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) << "MODE_TWO_FRACTURE_ENERGY is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS) || rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR)) << "TENSILE_INTERFACE_MODULUS is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS) || rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR)) << "SHEAR_INTERFACE_MODULUS is not a defined value" << std::endl;

    // Check the size of the vectors
    if (rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "INTERFACIAL_NORMAL_STRENGTH_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "INTERFACIAL_SHEAR_STRENGTH_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "MODE_ONE_FRACTURE_ENERGY_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "MODE_TWO_FRACTURE_ENERGY_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[TENSILE_INTERFACE_MODULUS_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "TENSILE_INTERFACE_MODULUS_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_INTERFACE_MODULUS_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "SHEAR_INTERFACE_MODULUS_VECTOR badly defined" << std::endl;
    }

    // Check negative fracture energy
    if (rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) {
        const double SizeModeOne = rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR].size();
        for (IndexType i=0; i < SizeModeOne; ++i) {
            KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] < 0.0) << "MODE_ONE_FRACTURE_ENERGY is negative at interface " << i << std::endl;
        }
    } else {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY] < 0.0) << "MODE_ONE_FRACTURE_ENERGY is negative." << std::endl;
    }

    if (rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) {
        const double SizeModeTwo = rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR].size();
        for (IndexType i=0; i < SizeModeTwo; ++i) {
            KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] < 0.0) << "MODE_TWO_FRACTURE_ENERGY is negative at interface " << i << std::endl;
        }
    } else {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY] < 0.0) << "MODE_TWO_FRACTURE_ENERGY is negative." << std::endl;
    }

    // Check fracture energy
    const double characteristic_length = 0.6343 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rElementGeometry));

    for(IndexType i=0; i < (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1; ++i) {
        const double Tn = rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
        const double GI = rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : rMaterialProperties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
        const double E = rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? rMaterialProperties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : rMaterialProperties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface

        const double AParameter_mode_one = 1.00 / (GI * E / (characteristic_length * std::pow(Tn, 2)) - 0.5); // Exponential
        KRATOS_ERROR_IF(AParameter_mode_one < 0.0) << "MODE_ONE_FRACTURE_ENERGY is too low at interface " << i << std::endl;
    }

    for(IndexType i=0; i < (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1; ++i) {
        const double Ts = rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
        const double GII = rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : rMaterialProperties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
        const double G = rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? rMaterialProperties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : rMaterialProperties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

        const double AParameter_mode_two = 1.00 / (GII * G / (characteristic_length * std::pow(Ts, 2)) - 0.5); // Exponential
        KRATOS_ERROR_IF(AParameter_mode_two < 0.0) << "MODE_TWO_FRACTURE_ENERGY is too low at interface " << i << std::endl;
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 4);
    } else {
        BaseType::CalculateTangentTensor(rValues,rStressMeasure);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class TractionSeparationLaw3D<3>;

} // Namespace Kratos
