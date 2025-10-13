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

ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::ThicknessIntegratedIsotropicConstitutiveLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::ThicknessIntegratedIsotropicConstitutiveLaw(
    const IndexType& rThicknessIntegrationPoints
    ) :
    ConstitutiveLaw()
{
    KRATOS_ERROR_IF(rThicknessIntegrationPoints <= 0) << "Wrong number of integration points through the thickness... " << std::endl;

    if (rThicknessIntegrationPoints != 5) // 5 is the default
        mThicknessIntegrationPoints = rThicknessIntegrationPoints;
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::ThicknessIntegratedIsotropicConstitutiveLaw(
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

ConstitutiveLaw::Pointer ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Create(
    Kratos::Parameters NewParameters
) const
{

    if (NewParameters.Has("thickness_integration_points")) {
        const IndexType thickness_integration_points = NewParameters["thickness_integration_points"].GetInt();
        KRATOS_ERROR_IF(thickness_integration_points <= 0) << "Wrong number of integration points through the thickness... " << std::endl;
        mThicknessIntegrationPoints = thickness_integration_points;
    }

    // We create the law
    return Kratos::make_shared<ThicknessIntegratedIsotropicConstitutiveLaw>(combination_factors);
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
    IndexType counter = 0;
    SizeType dimension = 3;
    if (mConstitutiveLaws.size() == 0) return dimension; // In case of not initialized CL
    // We perform the check in each layer
    for (auto& p_law : mConstitutiveLaws) {
        if (counter == 0) {
            dimension = p_law->WorkingSpaceDimension();
        } else {
            KRATOS_ERROR_IF_NOT(dimension == p_law->WorkingSpaceDimension()) << "Combining different size laws" << std::endl;
        }

        ++counter;
    }

    return dimension;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedIsotropicConstitutiveLaw::GetStrainSize() const
{
    IndexType counter = 0;
    SizeType strain_size = 6;
    if (mConstitutiveLaws.size() == 0) return strain_size; // In case of not initialized CL
    // We perform the check in each layer
    for (auto& p_law : mConstitutiveLaws) {
        if (counter == 0) {
            strain_size = p_law->GetStrainSize();
        } else {
            KRATOS_ERROR_IF_NOT(strain_size == p_law->GetStrainSize()) << "Combining different size laws" << std::endl;
        }

        ++counter;
    }

    return strain_size;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<bool>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<int>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<double>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<Matrix>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    // At least one layer should have the value
    rValue = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->GetValue(rThisVariable, rValue))
            break;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    // At least one layer should have the value
    rValue = 0;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            p_law->GetValue(rThisVariable, rValue);
            break;
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    rValue = 0.0;
    double aux_value;
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        const double factor = mCombinationFactors[i_layer];

        // we average over the layers
        if (p_law->Has(rThisVariable)) {
            p_law->GetValue(rThisVariable, aux_value);
            rValue += aux_value * factor;
        }
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    rValue.resize(VoigtSize, false);
    rValue.clear();
    Vector aux_value;
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        if (p_law->Has(rThisVariable)) {
            p_law->GetValue(rThisVariable, aux_value);
            rValue += aux_value * factor;
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Matrix& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        Matrix aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
array_1d<double, 3 >& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(3);
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        array_1d<double, 3 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
array_1d<double, 6 >& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(6);
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        array_1d<double, 6 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers

    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& p_law : mConstitutiveLaws) {
        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();
    const Vector strain_vector = rParameterValues.GetStrainVector();
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

    // We combine the value of each layer
    rValue = 0.0;
    double aux_value = 0.0;
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        Properties& r_prop = *(it_prop_begin + i_layer);
        rParameterValues.SetMaterialProperties(r_prop);

        // We rotate to local axes the strain
        noalias(rParameterValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

        aux_value = 0.0;
        p_law->CalculateValue(rParameterValues, rThisVariable, aux_value);
        rValue += factor * aux_value;

        // We reset the strain to its original global axes
        noalias(rParameterValues.GetStrainVector()) = strain_vector;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(r_material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We do some special operation for strains and stresses
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

        // We compute the strain
        if (rThisVariable == STRAIN) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == ALMANSI_STRAIN_VECTOR) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        if (rThisVariable == STRESSES) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        } if (rThisVariable == CAUCHY_STRESS_VECTOR) {
            this->CalculateMaterialResponseCauchy(rParameterValues);
        } if (rThisVariable == PK2_STRESS_VECTOR) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else {
        if (rThisVariable != DELAMINATION_DAMAGE_VECTOR_MODE_ONE && rThisVariable != DELAMINATION_DAMAGE_VECTOR_MODE_TWO) {
        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();
        const Vector strain_vector = rParameterValues.GetStrainVector();

        // We combine the value of each layer
        rValue.resize(VoigtSize, false);
        rValue.clear();
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
        for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
            const double factor = mCombinationFactors[i_layer];
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
            Properties& r_prop = *(it_prop_begin + i_layer);

            // We rotate to local axes the strain
            noalias(rParameterValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rParameterValues.SetMaterialProperties(r_prop);
            Vector aux_value;
            p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);

            // we return the aux_value to the global coordinates
            aux_value = prod(trans(voigt_rotation_matrix), aux_value);

            noalias(rValue) += factor * aux_value;

            noalias(rParameterValues.GetStrainVector()) = strain_vector;
        }

        // Reset properties
        rParameterValues.SetMaterialProperties(r_material_properties);
        }
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Matrix& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateValue(
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
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, false );

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
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else {
        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

        // We combine the value of each layer
        rValue.clear();
        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
        for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
            const double factor = mCombinationFactors[i_layer];
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
            Properties& r_prop = *(it_prop_begin + i_layer);

            rParameterValues.SetMaterialProperties(r_prop);
            Matrix aux_value;
            p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
            noalias(rValue) += factor * aux_value;
        }

        // Reset properties
        rParameterValues.SetMaterialProperties(r_material_properties);
    }

    return( rValue );
}
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
array_1d<double, 3 >& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    noalias(rValue) = ZeroVector(3);
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        Properties& r_prop = *(it_prop_begin + i_layer);

        rParameterValues.SetMaterialProperties(r_prop);
        array_1d<double, 3 > aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        noalias(rValue) += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(r_material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
array_1d<double, 6 >& ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    noalias(rValue) = ZeroVector(6);
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        Properties& r_prop = *(it_prop_begin + i_layer);

        rParameterValues.SetMaterialProperties(r_prop);
        array_1d<double, 6 > aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        noalias(rValue) += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(r_material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it layer by layer
    bool valid_input = true;
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
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
ConstitutiveLaw::StrainMeasure ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetStrainMeasure()
{
    // We return the first one
    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No constitutive laws defined" << std::endl;
    return mConstitutiveLaws[0]->GetStrainMeasure();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::StressMeasure ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetStressMeasure()
{
    // We return the first one
    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No constitutive laws defined" << std::endl;
    return mConstitutiveLaws[0]->GetStressMeasure();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::IsIncremental()
{
    // We check it layer by layer
    bool is_incremental = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->IsIncremental()) {
            is_incremental = true;
            break;
        }
    }

    return is_incremental;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // Resizing first
    mConstitutiveLaws.resize(mCombinationFactors.size());

    // We create the inner constitutive laws
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        Properties& r_prop = *(it_cl_begin + i_layer);

        KRATOS_ERROR_IF_NOT(r_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
        mConstitutiveLaws[i_layer] = r_prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLaws[i_layer]->InitializeMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    }

    KRATOS_DEBUG_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No CL defined" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);

    if (rValues.IsSetDeterminantF()) {
        Vector& r_stress_vector                = rValues.GetStressVector();
        const Matrix& r_deformation_gradient_f = rValues.GetDeformationGradientF();
        const double determinant_f             = rValues.GetDeterminantF();

        TransformStresses(r_stress_vector, r_deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Auxiliary stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();
        Vector auxiliary_stress_vector  = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {

            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
            const double factor            = mCombinationFactors[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(auxiliary_stress_vector) += factor * rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }
        noalias(rValues.GetStressVector()) = auxiliary_stress_vector;

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2);
        }

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Auxiliary stress vector
        const auto it_prop_begin            = r_material_properties.GetSubProperties().begin();
        Vector auxiliary_stress_vector       = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {

            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
            const double factor            = mCombinationFactors[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(auxiliary_stress_vector) += factor * rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }
        Vector &r_stress_vector = rValues.GetStressVector();
        noalias(r_stress_vector) = auxiliary_stress_vector;

        if (rValues.IsSetDeterminantF()) {
            // we push forward the stress
            Matrix stress_matrix(Dimension, Dimension);
            noalias(stress_matrix) = MathUtils<double>::StressVectorToTensor(r_stress_vector);
            ContraVariantPushForward(stress_matrix, rValues.GetDeformationGradientF()); // Kirchhoff
            noalias(r_stress_vector) = MathUtils<double>::StressTensorToVector( stress_matrix, r_stress_vector.size() );
        }

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2);
            // push forward Constitutive tangent tensor
            if (rValues.IsSetDeterminantF())
                PushForwardConstitutiveMatrix(rValues.GetConstitutiveMatrix(), rValues.GetDeformationGradientF());
        }

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);

    // Set to Cauchy Stress:
    if (rValues.IsSetDeterminantF()) {
        Vector& r_stress_vector       = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        const double determinant_f    = rValues.GetDeterminantF();

        r_stress_vector       /= determinant_f;
        r_constitutive_matrix /= determinant_f;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::InitializeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->InitializeMaterialResponsePK1(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::InitializeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->InitializeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->InitializeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->InitializeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->FinalizeMaterialResponsePK1(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->FinalizeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->FinalizeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }
    // The rotation matrix
    BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    const Vector strain_vector = rValues.GetStrainVector();
    // We perform the reset in each layer
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
        Properties& r_prop             = *(it_prop_begin + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        rValues.SetMaterialProperties(r_prop);
        // We rotate to local axes the strain
        noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
        p_law->FinalizeMaterialResponsePK2(rValues);
    }
    rValues.SetMaterialProperties(r_material_properties);
    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We perform the reset in each layer
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->ResetMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    }
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;

    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No constitutive laws defined" << std::endl;

    // We perform the check in each layer
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        aux_out += p_law->Check(r_prop, rElementGeometry, rCurrentProcessInfo);
    }

    if (rMaterialProperties.Has(LAYER_EULER_ANGLES)) {
        KRATOS_ERROR_IF(rMaterialProperties[LAYER_EULER_ANGLES].size() != 3 * mConstitutiveLaws.size()) << "Euler angles badly defined" << std::endl;
    }

    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateRotationMatrix(
        const Properties& rMaterialProperties,
        BoundedMatrix<double, VoigtSize, VoigtSize>& rRotationMatrix,
        const IndexType Layer
    )
{

    if (rRotationMatrix.size1() != VoigtSize)
        rRotationMatrix.resize(VoigtSize, VoigtSize, false);

    if (rMaterialProperties.Has(LAYER_EULER_ANGLES)) {
        const Vector layers_euler_angles = rMaterialProperties[LAYER_EULER_ANGLES];
        const double euler_angle_phi     = layers_euler_angles[3*Layer];
        const double euler_angle_theta   = layers_euler_angles[3*Layer + 1];
        const double euler_angle_hi      = layers_euler_angles[3*Layer + 2];

        BoundedMatrix<double, Dimension, Dimension>  rotation_matrix;

        if (std::abs(euler_angle_phi) + std::abs(euler_angle_theta) + std::abs(euler_angle_hi) > machine_tolerance) {
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(euler_angle_phi, euler_angle_theta, euler_angle_hi, rotation_matrix);
            ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(rotation_matrix, rRotationMatrix);
        } else {
            noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
        }
    } else {
        noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateAlmansiStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliary values
    const SizeType dimension = WorkingSpaceDimension();
    Vector& r_strain_vector = rValues.GetStrainVector();

    Matrix F(dimension, dimension);
    noalias(F) = rValues.GetDeformationGradientF();
    Matrix B_tensor;
    B_tensor.resize(dimension, dimension, false);
    noalias(B_tensor) = prod(F, trans(F));

    AdvancedConstitutiveLawUtilities<6>::CalculateAlmansiStrain(B_tensor, r_strain_vector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void ThicknessIntegratedIsotropicConstitutiveLaw<TDim>::CalculateGreenLagrangeStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliary values
    const SizeType dimension = WorkingSpaceDimension();
    Vector& r_strain_vector = rValues.GetStrainVector();

    Matrix F(dimension, dimension);
    noalias(F) = rValues.GetDeformationGradientF();
    Matrix C_tensor;
    C_tensor.resize(dimension, dimension, false);
    noalias(C_tensor) = prod(trans(F),F);

    ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, r_strain_vector);
}

/***********************************************************************************/
/***********************************************************************************/

template class ThicknessIntegratedIsotropicConstitutiveLaw<2>;
template class ThicknessIntegratedIsotropicConstitutiveLaw<3>;

} // Namespace Kratos
