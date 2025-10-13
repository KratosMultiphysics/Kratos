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
//  Main authors:    Alejandro Cornejo
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "thickness_integrated_isotropic_constitutive_law.h"
// #include "constitutive_laws_application_variables.h"
// #include "custom_utilities/tangent_operator_calculator_utility.h"
// #include "custom_utilities/advanced_constitutive_law_utilities.h"
// #include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw(
    const IndexType rThicknessIntegrationPoints
    ) :
    ConstitutiveLaw()
{
    KRATOS_ERROR_IF(rThicknessIntegrationPoints <= 0) << "Wrong number of integration points through the thickness... " << std::endl;

    if (rThicknessIntegrationPoints != 5) // 5 is the default
        mThicknessIntegrationPoints = rThicknessIntegrationPoints;
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw(
    const ThicknessIntegratedIsotropicConstitutiveLaw &rOther)
    : ConstitutiveLaw(rOther),
      mConstitutiveLaws(rOther.mConstitutiveLaws),
      mThicknessIntegrationPoints(rOther.mThicknessIntegrationPoints)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedIsotropicConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<ThicknessIntegratedIsotropicConstitutiveLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedIsotropicConstitutiveLaw::Create(
    Kratos::Parameters NewParameters
) const
{
    IndexType thickness_integration_points = 5; // Default value
    // We check if the user has defined a different value
    if (NewParameters.Has("thickness_integration_points")) {
        thickness_integration_points = NewParameters["thickness_integration_points"].GetInt();
        KRATOS_ERROR_IF(thickness_integration_points <= 0) << "Wrong number of integration points through the thickness... " << std::endl;
    }

    // We create the law
    return Kratos::make_shared<ThicknessIntegratedIsotropicConstitutiveLaw>(thickness_integration_points);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::~ThicknessIntegratedIsotropicConstitutiveLaw()
{
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedIsotropicConstitutiveLaw::WorkingSpaceDimension()
{

    return Dimension; // 3
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedIsotropicConstitutiveLaw::GetStrainSize() const
{
    return VoigtSize; // 8
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<bool>& rThisVariable)
{
    return Has<bool>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<int>& rThisVariable)
{
    return Has<int>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<double>& rThisVariable)
{
    return Has<double>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
{
    return Has<Vector>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return Has<Matrix>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    // for (auto& p_law : mConstitutiveLaws) {
    //     if (p_law->Has(rThisVariable)) {
    //         has = true;
    //         break;
    //     }
    // }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    // for (auto& p_law : mConstitutiveLaws) {
    //     if (p_law->Has(rThisVariable)) {
    //         has = true;
    //         break;
    //     }
    // }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    // At least one layer should have the value
    rValue = false;

    // for (auto& p_law : mConstitutiveLaws) {
    //     if (p_law->GetValue(rThisVariable, rValue))
    //         break;
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    // At least one layer should have the value
    rValue = 0;

    // for (auto& p_law : mConstitutiveLaws) {
    //     if (p_law->Has(rThisVariable)) {
    //         p_law->GetValue(rThisVariable, rValue);
    //         break;
    //     }
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    // rValue = 0.0;
    // double aux_value;
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     const double factor = mCombinationFactors[i_layer];

    //     // we average over the layers
    //     if (p_law->Has(rThisVariable)) {
    //         p_law->GetValue(rThisVariable, aux_value);
    //         rValue += aux_value * factor;
    //     }
    // }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    // rValue.resize(VoigtSize, false);
    // rValue.clear();
    // Vector aux_value;
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     if (p_law->Has(rThisVariable)) {
    //         p_law->GetValue(rThisVariable, aux_value);
    //         rValue += aux_value * factor;
    //     }
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    // rValue.clear();
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     Matrix aux_value;
    //     p_law->GetValue(rThisVariable, aux_value);
    //     rValue += aux_value * factor;
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    // rValue = ZeroVector(3);
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     array_1d<double, 3 > aux_value;
    //     p_law->GetValue(rThisVariable, aux_value);
    //     rValue += aux_value * factor;
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    // rValue = ZeroVector(6);
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     array_1d<double, 6 > aux_value;
    //     p_law->GetValue(rThisVariable, aux_value);
    //     rValue += aux_value * factor;
    // }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers

    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    // for (auto& p_law : mConstitutiveLaws) {
    //     p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    // }
}


/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();
    // const Vector strain_vector = rParameterValues.GetStrainVector();
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

    // // We combine the value of each layer
    // rValue = 0.0;
    // double aux_value = 0.0;
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     Properties& r_prop = *(it_prop_begin + i_layer);
    //     rParameterValues.SetMaterialProperties(r_prop);

    //     // We rotate to local axes the strain
    //     noalias(rParameterValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

    //     aux_value = 0.0;
    //     p_law->CalculateValue(rParameterValues, rThisVariable, aux_value);
    //     rValue += factor * aux_value;

    //     // We reset the strain to its original global axes
    //     noalias(rParameterValues.GetStrainVector()) = strain_vector;
    // }

    // // Reset properties
    // rParameterValues.SetMaterialProperties(r_material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return (rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return (rValue);
}
/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // // Resizing first
    // mConstitutiveLaws.resize(mCombinationFactors.size());

    // // We create the inner constitutive laws
    // const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     Properties& r_prop = *(it_cl_begin + i_layer);

    //     KRATOS_ERROR_IF_NOT(r_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    //     mConstitutiveLaws[i_layer] = r_prop[CONSTITUTIVE_LAW]->Clone();
    //     mConstitutiveLaws[i_layer]->InitializeMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    // }

    // KRATOS_DEBUG_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No CL defined" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{

}

/***********************************************************************************/
/***********************************************************************************/


void  ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{

}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK1(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK1(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}


/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We perform the reset in each layer
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     p_law->ResetMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    // }
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 8; // 3 membrane, 3 bending, 2 shear

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/


int ThicknessIntegratedIsotropicConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;

    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No constitutive laws defined" << std::endl;

    // We perform the check in each layer
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     aux_out += p_law->Check(r_prop, rElementGeometry, rCurrentProcessInfo);
    // }

    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/


// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateRotationMatrix(
//         const Properties& rMaterialProperties,
//         BoundedMatrix<double, VoigtSize, VoigtSize>& rRotationMatrix,
//         const IndexType Layer
//     )
// {
    // if (rRotationMatrix.size1() != VoigtSize)
    //     rRotationMatrix.resize(VoigtSize, VoigtSize, false);

    // if (rMaterialProperties.Has(LAYER_EULER_ANGLES)) {
    //     const Vector layers_euler_angles = rMaterialProperties[LAYER_EULER_ANGLES];
    //     const double euler_angle_phi     = layers_euler_angles[3*Layer];
    //     const double euler_angle_theta   = layers_euler_angles[3*Layer + 1];
    //     const double euler_angle_hi      = layers_euler_angles[3*Layer + 2];

    //     BoundedMatrix<double, Dimension, Dimension>  rotation_matrix;

    //     if (std::abs(euler_angle_phi) + std::abs(euler_angle_theta) + std::abs(euler_angle_hi) > machine_tolerance) {
    //         AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(euler_angle_phi, euler_angle_theta, euler_angle_hi, rotation_matrix);
    //         ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(rotation_matrix, rRotationMatrix);
    //     } else {
    //         noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
    //     }
    // } else {
    //     noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
    // }
// }

/***********************************************************************************/
/***********************************************************************************/

// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateAlmansiStrain(ConstitutiveLaw::Parameters& rValues)
// {
//     // Some auxiliary values
//     const SizeType dimension = WorkingSpaceDimension();
//     Vector& r_strain_vector = rValues.GetStrainVector();

//     Matrix F(dimension, dimension);
//     noalias(F) = rValues.GetDeformationGradientF();
//     Matrix B_tensor;
//     B_tensor.resize(dimension, dimension, false);
//     noalias(B_tensor) = prod(F, trans(F));

//     AdvancedConstitutiveLawUtilities<6>::CalculateAlmansiStrain(B_tensor, r_strain_vector);
// }

/***********************************************************************************/
/***********************************************************************************/


// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateGreenLagrangeStrain(ConstitutiveLaw::Parameters& rValues)
// {
//     // Some auxiliary values
//     const SizeType dimension = WorkingSpaceDimension();
//     Vector& r_strain_vector = rValues.GetStrainVector();

//     Matrix F(dimension, dimension);
//     noalias(F) = rValues.GetDeformationGradientF();
//     Matrix C_tensor;
//     C_tensor.resize(dimension, dimension, false);
//     noalias(C_tensor) = prod(trans(F),F);

//     ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, r_strain_vector);
// }

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
