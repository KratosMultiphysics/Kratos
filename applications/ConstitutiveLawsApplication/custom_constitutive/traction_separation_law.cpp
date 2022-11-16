// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Fernando Rastellini
//                   Alejandro Cornejo Velazquez
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/traction_separation_law.h"
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
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const std::vector<double>& rCombinationFactors) : BaseType()
{
    // We compute the proportion of the factors (must be over 1)
    double aux_factor = 0.0;
    for (IndexType i_layer = 0; i_layer < rCombinationFactors.size(); ++i_layer) {
        aux_factor += rCombinationFactors[i_layer];
    }

    KRATOS_ERROR_IF(aux_factor < std::numeric_limits<double>::epsilon()) << "Wrong factors in TractionSeparationLaw3D" << std::endl;

    // Resize
    mCombinationFactors.resize(rCombinationFactors.size());

    // We fill the maps
    for (IndexType i_layer = 0; i_layer < rCombinationFactors.size(); ++i_layer) {
        mCombinationFactors[i_layer] = rCombinationFactors[i_layer]/aux_factor;
    }
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const TractionSeparationLaw3D<TDim>& rOther)
    : BaseType(rOther),
      mConstitutiveLaws(rOther.mConstitutiveLaws),
      mCombinationFactors(rOther.mCombinationFactors),
    //   mGc(rOther.mGc),
    //   minitial_threshold(rOther.minitial_threshold),
    //   mthreshold(rOther.mthreshold),
      mdelamination_damage(rOther.mdelamination_damage),
    //   mAParameter(rOther.mAParameter),
    //   mDamageIndicator(rOther.mDamageIndicator),
      mstatus_coeff(rOther.mstatus_coeff)
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
std::size_t TractionSeparationLaw3D<TDim>::WorkingSpaceDimension()
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

template<unsigned int TDim>
std::size_t TractionSeparationLaw3D<TDim>::GetStrainSize() const
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<bool>& rThisVariable)
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<int>& rThisVariable)
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<double>& rThisVariable)
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;

    for (auto& p_law : mConstitutiveLaws) {
        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR) {
        return true;
    }
    
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<Matrix>& rThisVariable)
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
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
bool TractionSeparationLaw3D<TDim>::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
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
bool& TractionSeparationLaw3D<TDim>::GetValue(
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
int& TractionSeparationLaw3D<TDim>::GetValue(
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
double& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    rValue = 0.0;
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        double aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    // for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
    //     const double factor = mCombinationFactors[i_layer];
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     Vector aux_value;
    //     p_law->GetValue(rThisVariable, aux_value);
    //     rValue += aux_value * factor;
    // }

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR) {
        
        rValue.resize(mCombinationFactors.size()+1, false);
        
        noalias(rValue) = mdelamination_damage;
        return rValue;
    }
    
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Matrix& TractionSeparationLaw3D<TDim>::GetValue(
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
array_1d<double, 3 >& TractionSeparationLaw3D<TDim>::GetValue(
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
array_1d<double, 6 >& TractionSeparationLaw3D<TDim>::GetValue(
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
void TractionSeparationLaw3D<TDim>::SetValue(
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
void TractionSeparationLaw3D<TDim>::SetValue(
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
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& TractionSeparationLaw3D<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = 0.0;
    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mCombinationFactors.size(); ++i_layer) {
        const double factor = mCombinationFactors[i_layer];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        Properties& r_prop = *(it_prop_begin + i_layer);

        rParameterValues.SetMaterialProperties(r_prop);
        double aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(r_material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& TractionSeparationLaw3D<TDim>::CalculateValue(
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
        const bool flag_const_tensor = r_flags.Is( BaseType::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( BaseType::COMPUTE_STRESS );

        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( BaseType::COMPUTE_STRESS, false );

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
        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( BaseType::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( BaseType::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( BaseType::COMPUTE_STRESS );

        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( BaseType::COMPUTE_STRESS, true );

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
        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( BaseType::COMPUTE_STRESS, flag_stress );
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
            Vector aux_value;
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
Matrix& TractionSeparationLaw3D<TDim>::CalculateValue(
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
        const bool flag_const_tensor = r_flags.Is( BaseType::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( BaseType::COMPUTE_STRESS );

        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( BaseType::COMPUTE_STRESS, false );

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
        r_flags.Set( BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( BaseType::COMPUTE_STRESS, flag_stress );
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
array_1d<double, 3 >& TractionSeparationLaw3D<TDim>::CalculateValue(
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
array_1d<double, 6 >& TractionSeparationLaw3D<TDim>::CalculateValue(
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
bool TractionSeparationLaw3D<TDim>::ValidateInput(const Properties& rMaterialProperties)
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
ConstitutiveLaw::StrainMeasure TractionSeparationLaw3D<TDim>::GetStrainMeasure()
{
    // We return the first one
    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "TractionSeparationLaw3D: No constitutive laws defined" << std::endl;
    return mConstitutiveLaws[0]->GetStrainMeasure();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::StressMeasure TractionSeparationLaw3D<TDim>::GetStressMeasure()
{
    // We return the first one
    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "TractionSeparationLaw3D: No constitutive laws defined" << std::endl;
    return mConstitutiveLaws[0]->GetStressMeasure();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::IsIncremental()
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
void TractionSeparationLaw3D<TDim>::InitializeMaterial(
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

    KRATOS_DEBUG_ERROR_IF(mConstitutiveLaws.size() == 0) << "TractionSeparationLaw3D: No CL defined" << std::endl;

    // mthreshold.resize(mConstitutiveLaws.size()-1, false);
    // for (int i=0; i < mConstitutiveLaws.size()-1; ++i) {
    //         mthreshold[i] = std::numeric_limits<double>::infinity();
    //     }
    
    // mDamageIndicator.resize(mConstitutiveLaws.size()-1, false);
    // for (int i=0; i < mConstitutiveLaws.size()-1; ++i) {
    //         mDamageIndicator[i] = 1;
    //     }

    mstatus_coeff.resize(mConstitutiveLaws.size()-1, false);
    for (int i=0; i < mConstitutiveLaws.size()-1; ++i) {
            mstatus_coeff[i] = 1;
        }
    
    // mGc.resize(mConstitutiveLaws.size()-1, false);
    // noalias(mGc) = ZeroVector(mConstitutiveLaws.size()-1);

    // minitial_threshold.resize(mConstitutiveLaws.size()-1, false);
    // noalias(minitial_threshold) = ZeroVector(mConstitutiveLaws.size()-1);

    mdelamination_damage.resize(mConstitutiveLaws.size()+1, false);
    noalias(mdelamination_damage) = ZeroVector(mConstitutiveLaws.size()+1);

    // mAParameter.resize(mConstitutiveLaws.size()-1, false);
    // noalias(mAParameter) = ZeroVector(mConstitutiveLaws.size()-1);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
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
        CalculateGreenLagrangeStrain(rValues);
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
        // Vector undamaged_auxiliar_stress_vector  = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;


        std::vector<Vector> layer_stress(mConstitutiveLaws.size());
        for (int i=0; i < mConstitutiveLaws.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> delamination_damage_affected_stress_matrix(mConstitutiveLaws.size());
        for (int i=0; i < mConstitutiveLaws.size(); ++i) {
            delamination_damage_affected_stress_matrix[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(mConstitutiveLaws.size()-1);
        for (int i=0; i < mConstitutiveLaws.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }


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
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();
            noalias(delamination_damage_affected_stress_matrix[i_layer]) = rValues.GetStressVector();
            // noalias(undamaged_auxiliar_stress_vector) += factor * rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        // Vector Gc(mConstitutiveLaws.size()-1);
        // Vector initial_threshold(mConstitutiveLaws.size()-1);
        // Vector threshold(mConstitutiveLaws.size()-1);
        Vector delamination_damage(mConstitutiveLaws.size()+1);
        Vector negative_interfacial_stress_index(mConstitutiveLaws.size()+1);
        noalias(negative_interfacial_stress_index) = ZeroVector(mConstitutiveLaws.size()+1);
        // Vector AParameter(mConstitutiveLaws.size()-1);
        // Vector DamageIndicator(mConstitutiveLaws.size()-1);
        Vector status_coeff(mConstitutiveLaws.size()-1);

        // noalias(Gc) = mGc;
        // noalias(initial_threshold) = minitial_threshold;
        // noalias(threshold) = mthreshold;
        noalias(delamination_damage) = mdelamination_damage;
        // noalias(AParameter) = mAParameter;
        // noalias(DamageIndicator) = mDamageIndicator;
        noalias(status_coeff) = mstatus_coeff;


        for(int i=0; i < mConstitutiveLaws.size()-1; ++i) {

            interfacial_stress[i][0] = (MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5)); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress
            
            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < 0.0) {
                negative_interfacial_stress_index [i+1] = 1.0;
            }

            // Von mises equivalent stress
            Vector interfacial_stress_vector = ZeroVector(VoigtSize);

            interfacial_stress_vector[2] = interfacial_stress[i][0];
            interfacial_stress_vector[4] = interfacial_stress[i][1];
            interfacial_stress_vector[5] = interfacial_stress[i][2];
            double I1, J2;
            array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);

            ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(interfacial_stress_vector, I1);
            ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(interfacial_stress_vector, I1, deviator, J2);
            // End Von mises equivalent stress

            // Damage calculation

            const double T0n = r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double T0t = r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Eta = r_material_properties[BK_COEFFICIENT]; // Benzeggagh-Kenane (B-K) Law Coefficient
            const double Ei = r_material_properties[TENSILE_INTERAFCE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties[SHEAR_INTERAFCE_MODULUS]; // Shear modulus of the interface
            // const double characteristic_length = 0.0003; // Characteristic Length of the Cohesive Part
            // const double characteristic_length = 0.25 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry()));
            const double characteristic_length = 0.000165;
            const double tolerance = std::numeric_limits<double>::epsilon();
            // const double Fd = std::pow(interfacial_stress[i][0]/T0n,2.0)+std::pow(interfacial_stress[i][1]/T0s,2.0)+std::pow(interfacial_stress[i][2]/T0t,2.0); // Damage Initiation Criterion
            double T_eq = std::sqrt(std::pow(interfacial_stress[i][0],2.0)+std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));
            // double T_eq = std::sqrt(3.0 * J2);
            const double T_shear = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));
            // const double mode_mix_factor = std::pow(T_shear / T_eq,2.0);
            const double beta_factor = T_shear / (T_shear + interfacial_stress[i][0]);
            const double mode_mix_factor = (beta_factor * beta_factor) / (1 + 2 * beta_factor * beta_factor - 2 * beta_factor);
            const double initial_threshold = std::sqrt(std::pow(T0n,2.0) + (std::pow(T0s,2.0) - std::pow(T0n,2.0)) * std::pow(mode_mix_factor, Eta));
            // const double initial_threshold = T0n;
            const double threshold = status_coeff[i] * initial_threshold;
            const double F = T_eq - threshold;
            if (F > tolerance) {
                const double Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
                // const double Gc = GIIc;
                const double K = Ei + (Gi - Ei) * std::pow(mode_mix_factor, Eta);
                // const double K = Gi;
                // const double AParameter = -std::pow(initial_threshold, 2) / (2.0 * K * Gc / characteristic_length); // Linear
                // const double AParameter = 1.00 / (Gc * K / (characteristic_length * std::pow(initial_threshold, 2)) - 0.5); // Exponential
                const double AParameter = 0.07;
                KRATOS_ERROR_IF(AParameter < 0.0) << "AParameter is negative." << std::endl;

                // delamination_damage[i+1] = (1.0 - initial_threshold / T_eq) / (1.0 + AParameter); // Linear
                delamination_damage[i+1] = 1.0 - (initial_threshold / T_eq) * std::exp(AParameter *
                    (1.0 - T_eq / initial_threshold)); // Exponential

                delamination_damage[i+1] = (delamination_damage[i+1] >= 0.99999) ? 0.99999 : delamination_damage[i+1];
                delamination_damage[i+1] = (delamination_damage[i+1] < 0.0) ? 0.0 : delamination_damage[i+1];

                // KRATOS_ERROR_IF(delamination_damage[i+1] < mdelamination_damage[i+1]) << "AParameter is negative." << std::endl;
            }

            // End damage calculation
        }

        for(int i=0; i < mConstitutiveLaws.size(); ++i) {
            double layer_damage_variable = 0;
            if (delamination_damage[i+1] > delamination_damage[i]) {
                layer_damage_variable = delamination_damage[i+1];
            } else {
                layer_damage_variable = delamination_damage[i];
            }
            delamination_damage_affected_stress_matrix[i][0] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][0];
            delamination_damage_affected_stress_matrix[i][1] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][1];
            delamination_damage_affected_stress_matrix[i][2] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][2];
            delamination_damage_affected_stress_matrix[i][3] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][3];
            delamination_damage_affected_stress_matrix[i][4] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][4];
            delamination_damage_affected_stress_matrix[i][5] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][5];

            // if (negative_interfacial_stress_index [i] != 1.0 && negative_interfacial_stress_index [i+1] != 1.0) {
            //     delamination_damage_affected_stress_matrix[i][2] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][2];
            // }
        }
        // Calculating output stresses

        // for(int i=0; i < mConstitutiveLaws.size(); ++i) {

        //     Properties& r_prop             = *(it_prop_begin + i);
        //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        //     const double factor            = mCombinationFactors[i];

        //     std::vector<double> layer_damage_vector(6);
        //     layer_damage_vector = {1,
        //                             1,
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1]),
        //                             1,
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1]),
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1])};

        //     auxiliar_stress_vector[0] += factor * layer_damage_vector[0] * layer_stress[i][0];
        //     auxiliar_stress_vector[1] += factor * layer_damage_vector[1] * layer_stress[i][1];
        //     auxiliar_stress_vector[2] += factor * layer_damage_vector[2] * layer_stress[i][2];
        //     auxiliar_stress_vector[3] += factor * layer_damage_vector[3] * layer_stress[i][3];
        //     auxiliar_stress_vector[4] += factor * layer_damage_vector[4] * layer_stress[i][4];
        //     auxiliar_stress_vector[5] += factor * layer_damage_vector[5] * layer_stress[i][5];
        // }

        // double damage_coeff = 1;

        for(int i=0; i < mConstitutiveLaws.size(); ++i) {
            const double factor = mCombinationFactors[i];
            delamination_damage_affected_stress_vector += factor * delamination_damage_affected_stress_matrix[i];
        }

        auxiliar_stress_vector = delamination_damage_affected_stress_vector;

        // auxiliar_stress_vector[0] = damage_coeff * undamaged_auxiliar_stress_vector[0];
        // auxiliar_stress_vector[1] = damage_coeff * undamaged_auxiliar_stress_vector[1];
        // auxiliar_stress_vector[2] = damage_coeff * undamaged_auxiliar_stress_vector[2];
        // auxiliar_stress_vector[3] = damage_coeff * undamaged_auxiliar_stress_vector[3];
        // auxiliar_stress_vector[4] = damage_coeff * undamaged_auxiliar_stress_vector[4];
        // auxiliar_stress_vector[5] = damage_coeff * undamaged_auxiliar_stress_vector[5];


        // End calculating output stresses

        //

        // Delamination Damage Criterion V1

        // const double T0n = 16000000; // Interfacial Normal Strength
        // const double T0s = 27000000; // Interfacial Shear Strength
        // const double T0t = 27000000; // Interfacial Shear Strength
        // const double GIc = 102; // Mode I Energy Release Rate
        // const double GIIc = 194; // Mode II Energy Release Rate
        // const double Eta = 2; // Benzeggagh-Kenane (B-K) Law Coefficient
        // double Gc = mGc; // Mix Mode Energy Release Rate
        // double Elastic_energy = mElastic_energy; // Elastic energy stored before damage initiation 
        // double Delta_G = mDelta_G;
        // double Delta_epsilon_one = mDelta_epsilon_one;
        // double SERR = mSERR; // Strain Energy Release Rate 
        // double delamination_damage = mdelamination_damage; // Scalar delamination damage variable  
        // double Delta_eq = 0; // Equivalent Strain
        // double Delta_eq_max = mDelta_eq_max; // Equivalent Strain History Variable
        // double T_eq = mT_eq; // Equivalent Stress
        // double DamageIndicator = mDamageIndicator; // Onset of Damage
        // // const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry()); // Characteristic Length of the Element
        // const double characteristic_length = 0.0001; // Characteristic Length of the Cohesive Part
        // double Fd = std::pow(undamaged_auxiliar_stress_vector[2]/T0n,2.0)+std::pow(undamaged_auxiliar_stress_vector[4]/T0s,2.0)+std::pow(undamaged_auxiliar_stress_vector[5]/T0t,2.0); // Damage Initiation Criterion
        // if (Fd >= 1.0 || DamageIndicator >= 2.0) {
        //     Delta_eq = std::sqrt(std::pow(strain_vector[2],2.0)+std::pow(strain_vector[4],2.0)+std::pow(strain_vector[5],2.0));
        //     if (DamageIndicator == 1.0) { // We calculate Elastic energy and mode mix only at damage initiation
        //         T_eq = std::sqrt(std::pow(undamaged_auxiliar_stress_vector[2],2.0)+std::pow(undamaged_auxiliar_stress_vector[4],2.0)+std::pow(undamaged_auxiliar_stress_vector[5],2.0));
        //         double Gn = undamaged_auxiliar_stress_vector[2] * strain_vector[2] / 2.0;
        //         double Gs = undamaged_auxiliar_stress_vector[4] * strain_vector[4] / 2.0;
        //         double Gt = undamaged_auxiliar_stress_vector[5] * strain_vector[5] / 2.0;
        //         double mode_mix_factor = (Gs+Gt) / (Gn+Gs+Gt);
        //         // Elastic_energy = (T_eq * Delta_eq /2.0);
        //         Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
        //         Delta_G = (Gc / characteristic_length) / (((2.0 * Gc / characteristic_length) / T_eq) - Delta_eq);
        //         // Delta_eq_max = Delta_eq;
        //         Delta_epsilon_one = Delta_eq;
        //     }
        //     if (Delta_eq >= Delta_eq_max) { // Loading
        //         // SERR += (Delta_eq - Delta_eq_max) * T_eq; // Strain Energy Release Rate
        //         // delamination_damage = std::min(SERR / ((Gc / characteristic_length) - Elastic_energy), 1.0);
        //         // Delta_eq_max = Delta_eq;

        //         SERR = (Delta_eq - Delta_epsilon_one) * Delta_G; // Strain Energy Release Rate
        //         delamination_damage = std::min(SERR / (Gc / characteristic_length),1.0);
        //         // Delta_eq_max = Delta_eq;
        //     }
        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5]; 
        //     // T_eq = std::sqrt(std::pow(auxiliar_stress_vector[2],2.0)+std::pow(auxiliar_stress_vector[4],2.0)+std::pow(auxiliar_stress_vector[5],2.0));
        //     // DamageIndicator += 1.0;
        // } else { // Undamaged Case
        //     auxiliar_stress_vector = undamaged_auxiliar_stress_vector;
        // }

        // End Delamination Damage Criterion V1

        //

        // Delamination Damage Criterion V2

        // const double T0n = 16000000; // Interfacial Normal Strength
        // const double T0s = 27000000; // Interfacial Shear Strength
        // const double T0t = 27000000; // Interfacial Shear Strength
        // const double GIc = 102; // Mode I Energy Release Rate
        // const double GIIc = 194; // Mode II Energy Release Rate
        // const double Eta = 2; // Benzeggagh-Kenane (B-K) Law Coefficient
        // double Gc = mGc; // Mix Mode Energy Release Rate
        // double initial_threshold = minitial_threshold;
        // double threshold = mthreshold;
        // double delamination_damage = mdelamination_damage;
        // double AParameter = mAParameter;
        // double DamageIndicator = mDamageIndicator;
        // const double characteristic_length = 0.0001; // Characteristic Length of the Cohesive Part
        // const double tolerance = std::numeric_limits<double>::epsilon();
        // const double Fd = std::pow(undamaged_auxiliar_stress_vector[2]/T0n,2.0)+std::pow(undamaged_auxiliar_stress_vector[4]/T0s,2.0)+std::pow(undamaged_auxiliar_stress_vector[5]/T0t,2.0); // Damage Initiation Criterion
        // double T_eq = std::sqrt(std::pow(undamaged_auxiliar_stress_vector[2],2.0)+std::pow(undamaged_auxiliar_stress_vector[4],2.0)+std::pow(undamaged_auxiliar_stress_vector[5],2.0));
        // if (Fd >= 1.0 && DamageIndicator == 1) {
        //     initial_threshold = T_eq;
        //     threshold = T_eq;
        //     double Gn = undamaged_auxiliar_stress_vector[2] * strain_vector[2] / 2.0;
        //     double Gs = undamaged_auxiliar_stress_vector[4] * strain_vector[4] / 2.0;
        //     double Gt = undamaged_auxiliar_stress_vector[5] * strain_vector[5] / 2.0;
        //     double mode_mix_factor = (Gs+Gt) / (Gn+Gs+Gt);
        //     Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
        //     // AParameter = -std::pow(initial_threshold, 2) / (2.0 * 800000000 * Gc / characteristic_length); // Linear
        //     AParameter = 1.00 / (Gc * 800000000 / (characteristic_length * std::pow(initial_threshold, 2)) - 0.5); // Exponential
        // }
        // const double Fp = T_eq - threshold;
        // if (Fp <= tolerance) { // Elastic case
        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5]; 
        // } else { // Damage case
        //     // delamination_damage = (1.0 - initial_threshold / T_eq) / (1.0 + AParameter); // Linear
        //     delamination_damage = 1.0 - (initial_threshold / T_eq) * std::exp(AParameter *
        //         (1.0 - T_eq / initial_threshold)); // Exponential

        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5];
        // }

        // End Delamination Damage Criterion V2

        //

        // Delamination Damage Criterion V3

        // const double T0n = 16000000; // Interfacial Normal Strength
        // const double T0s = 27000000; // Interfacial Shear Strength
        // const double T0t = 27000000; // Interfacial Shear Strength
        // const double GIc = 102; // Mode I Energy Release Rate
        // const double GIIc = 194; // Mode II Energy Release Rate
        // const double Eta = 2; // Benzeggagh-Kenane (B-K) Law Coefficient
        // double Gc = mGc; // Mix Mode Energy Release Rate
        // double initial_threshold = minitial_threshold;
        // double initial_Delta_eq = minitial_Delta_eq;
        // double threshold = 0;
        // double delamination_damage = mdelamination_damage;
        // double DamageIndicator = mDamageIndicator;
        // double Delta_eq_max = mDelta_eq_max;
        // const double characteristic_length = 0.0001; // Characteristic Length of the Cohesive Part
        // const double tolerance = std::numeric_limits<double>::epsilon();
        // double normal_stress = MacaullyBrackets(undamaged_auxiliar_stress_vector[2]);
        // double normal_strain = MacaullyBrackets(strain_vector[2]);
        // const double Fd = std::pow(normal_stress/T0n,2.0)+std::pow(undamaged_auxiliar_stress_vector[4]/T0s,2.0)+std::pow(undamaged_auxiliar_stress_vector[5]/T0t,2.0); // Damage Initiation Criterion
        // double T_eq = std::sqrt(std::pow(normal_stress,2.0)+std::pow(undamaged_auxiliar_stress_vector[4],2.0)+std::pow(undamaged_auxiliar_stress_vector[5],2.0));
        // double Delta_eq = std::sqrt(std::pow(normal_strain,2.0)+std::pow(strain_vector[4],2.0)+std::pow(strain_vector[5],2.0));
        // if (Fd >= 1.0 && DamageIndicator == 1) {
        //     initial_threshold = T_eq;
        //     initial_Delta_eq = Delta_eq;
        //     double Gn = normal_stress * normal_strain / 2.0;
        //     double Gs = undamaged_auxiliar_stress_vector[4] * strain_vector[4] / 2.0;
        //     double Gt = undamaged_auxiliar_stress_vector[5] * strain_vector[5] / 2.0;
        //     double mode_mix_factor = (Gs+Gt) / (Gn+Gs+Gt);
        //     Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
        // }
        // if (Delta_eq >= Delta_eq_max && DamageIndicator == 2) {
        //     threshold = initial_threshold * std::exp((initial_threshold*(Delta_eq-initial_Delta_eq))/(0.5*initial_threshold*initial_Delta_eq-(Gc / characteristic_length)));
        //     delamination_damage = 1 - (threshold / T_eq);
        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5];
        // } else {
        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5];
        // } 
        

        // End Delamination Damage Criterion V3

        //


        // KRATOS_WATCH(delamination_damage);
        
        noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, BaseType::StressMeasure_PK2);
        }

        // test for tangent tensor

        // const double E  = 800000000;
        // const double NU = 0.0;

        // const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        // const double c2 = c1 * ( 1 - NU );
        // const double c3 = c1 * NU;
        // const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

        // Matrix test_tangent_tensor = ZeroMatrix(6,6);

        // test_tangent_tensor( 0, 0 ) = c2;
        // test_tangent_tensor( 0, 1 ) = c3;
        // test_tangent_tensor( 0, 2 ) = c3;
        // test_tangent_tensor( 1, 0 ) = c3;
        // test_tangent_tensor( 1, 1 ) = c2;
        // test_tangent_tensor( 1, 2 ) = c3;
        // test_tangent_tensor( 2, 0 ) = c3;
        // test_tangent_tensor( 2, 1 ) = c3;
        // test_tangent_tensor( 2, 2 ) = c2;
        // test_tangent_tensor( 3, 3 ) = c4;
        // test_tangent_tensor( 4, 4 ) = c4;
        // test_tangent_tensor( 5, 5 ) = c4;

        // rValues.GetConstitutiveMatrix() = test_tangent_tensor;

        // test for tangent tensor

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
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is( BaseType::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin            = r_material_properties.GetSubProperties().begin();
        Vector auxiliar_stress_vector       = ZeroVector(VoigtSize);

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
            noalias(auxiliar_stress_vector) += factor * rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }
        Vector &r_stress_vector = rValues.GetStressVector();
        noalias(r_stress_vector) = auxiliar_stress_vector;

        if (rValues.IsSetDeterminantF()) {
            // we push forward the stress
            Matrix stress_matrix(Dimension, Dimension);
            noalias(stress_matrix) = MathUtils<double>::StressVectorToTensor(r_stress_vector);
            ContraVariantPushForward(stress_matrix, rValues.GetDeformationGradientF()); // Kirchhoff
            noalias(r_stress_vector) = MathUtils<double>::StressTensorToVector( stress_matrix, r_stress_vector.size() );
        }

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, BaseType::StressMeasure_PK2);
            // push forward Constitutive tangent tensor
            if (rValues.IsSetDeterminantF())
                PushForwardConstitutiveMatrix(rValues.GetConstitutiveMatrix(), rValues.GetDeformationGradientF());
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
void TractionSeparationLaw3D<TDim>::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
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
void TractionSeparationLaw3D<TDim>::InitializeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
void TractionSeparationLaw3D<TDim>::InitializeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
void TractionSeparationLaw3D<TDim>::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
void TractionSeparationLaw3D<TDim>::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
    r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
    r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
    // r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);


    //
    
    // Finalizing delamination damage


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
        CalculateGreenLagrangeStrain(rValues);
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
        // Vector undamaged_auxiliar_stress_vector  = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;


        std::vector<Vector> layer_stress(mConstitutiveLaws.size());
        for (int i=0; i < mConstitutiveLaws.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> delamination_damage_affected_stress_matrix(mConstitutiveLaws.size());
        for (int i=0; i < mConstitutiveLaws.size(); ++i) {
            delamination_damage_affected_stress_matrix[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(mConstitutiveLaws.size()-1);
        for (int i=0; i < mConstitutiveLaws.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }


        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {

            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
            const double factor            = mCombinationFactors[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);
            // p_law->FinalizeMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();
            noalias(delamination_damage_affected_stress_matrix[i_layer]) = rValues.GetStressVector();
            // noalias(undamaged_auxiliar_stress_vector) += factor * rValues.GetStressVector();

            p_law->FinalizeMaterialResponsePK2(rValues);

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        // Vector Gc(mConstitutiveLaws.size()-1);
        // Vector initial_threshold(mConstitutiveLaws.size()-1);
        // Vector threshold(mConstitutiveLaws.size()-1);
        Vector delamination_damage(mConstitutiveLaws.size()+1);
        Vector negative_interfacial_stress_index(mConstitutiveLaws.size()+1);
        noalias(negative_interfacial_stress_index) = ZeroVector(mConstitutiveLaws.size()+1);
        // Vector AParameter(mConstitutiveLaws.size()-1);
        // Vector DamageIndicator(mConstitutiveLaws.size()-1);
        Vector status_coeff(mConstitutiveLaws.size()-1);

        // noalias(Gc) = mGc;
        // noalias(initial_threshold) = minitial_threshold;
        // noalias(threshold) = mthreshold;
        noalias(delamination_damage) = mdelamination_damage;
        // noalias(AParameter) = mAParameter;
        // noalias(DamageIndicator) = mDamageIndicator;
        noalias(status_coeff) = mstatus_coeff;


        for(int i=0; i < mConstitutiveLaws.size()-1; ++i) {

            interfacial_stress[i][0] = (MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5)); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < 0.0) {
                negative_interfacial_stress_index [i+1] = 1.0;
            }

            // Von mises equivalent stress
            Vector interfacial_stress_vector = ZeroVector(VoigtSize);

            interfacial_stress_vector[2] = interfacial_stress[i][0];
            interfacial_stress_vector[4] = interfacial_stress[i][1];
            interfacial_stress_vector[5] = interfacial_stress[i][2];
            double I1, J2;
            array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);

            ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(interfacial_stress_vector, I1);
            ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(interfacial_stress_vector, I1, deviator, J2);
            // End Von mises equivalent stress

            // Damage calculation

            const double T0n = r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double T0t = r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Eta = r_material_properties[BK_COEFFICIENT]; // Benzeggagh-Kenane (B-K) Law Coefficient
            const double Ei = r_material_properties[TENSILE_INTERAFCE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties[SHEAR_INTERAFCE_MODULUS]; // Shear modulus of the interface
            // const double characteristic_length = 0.0003; // Characteristic Length of the Cohesive Part
            // const double characteristic_length = 0.25 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry()));
            const double characteristic_length = 0.000165;
            const double tolerance = std::numeric_limits<double>::epsilon();
            // const double Fd = std::pow(interfacial_stress[i][0]/T0n,2.0)+std::pow(interfacial_stress[i][1]/T0s,2.0)+std::pow(interfacial_stress[i][2]/T0t,2.0); // Damage Initiation Criterion
            double T_eq = std::sqrt(std::pow(interfacial_stress[i][0],2.0)+std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));
            // double T_eq = std::sqrt(3.0 * J2);
            const double T_shear = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));
            // const double mode_mix_factor = std::pow(T_shear / T_eq,2.0);
            const double beta_factor = T_shear / (T_shear + interfacial_stress[i][0]);
            const double mode_mix_factor = (beta_factor * beta_factor) / (1 + 2 * beta_factor * beta_factor - 2 * beta_factor);
            const double initial_threshold = std::sqrt(std::pow(T0n,2.0) + (std::pow(T0s,2.0) - std::pow(T0n,2.0)) * std::pow(mode_mix_factor, Eta));
            // const double initial_threshold = T0n;
            const double threshold = status_coeff[i] * initial_threshold;
            const double F = T_eq - threshold;
            if (F > tolerance) {
                const double Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
                // const double Gc = GIIc; // Benzeggagh-Kenane (B-K) Law
                const double K = Ei + (Gi - Ei) * std::pow(mode_mix_factor, Eta);
                // const double K = Gi;
                // const double AParameter = -std::pow(initial_threshold, 2) / (2.0 * K * Gc / characteristic_length); // Linear
                // const double AParameter = 1.00 / (Gc * K / (characteristic_length * std::pow(initial_threshold, 2)) - 0.5); // Exponential
                const double AParameter = 0.07;
                KRATOS_ERROR_IF(AParameter < 0.0) << "AParameter is negative." << std::endl;

                // delamination_damage[i+1] = (1.0 - initial_threshold / T_eq) / (1.0 + AParameter); // Linear
                delamination_damage[i+1] = 1.0 - (initial_threshold / T_eq) * std::exp(AParameter *
                    (1.0 - T_eq / initial_threshold)); // Exponential

                delamination_damage[i+1] = (delamination_damage[i+1] >= 0.99999) ? 0.99999 : delamination_damage[i+1]; //0.98
                delamination_damage[i+1] = (delamination_damage[i+1] < 0.0) ? 0.0 : delamination_damage[i+1];

                mdelamination_damage[i+1] = delamination_damage[i+1];
                mstatus_coeff[i] = T_eq / initial_threshold;
            }

            // End damage calculation
        }

        for(int i=0; i < mConstitutiveLaws.size(); ++i) {
            double layer_damage_variable = 0;
            if (delamination_damage[i+1] > delamination_damage[i]) {
                layer_damage_variable = delamination_damage[i+1];
            } else {
                layer_damage_variable = delamination_damage[i];
            }
            delamination_damage_affected_stress_matrix[i][0] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][0];
            delamination_damage_affected_stress_matrix[i][1] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][1];
            delamination_damage_affected_stress_matrix[i][2] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][2];
            delamination_damage_affected_stress_matrix[i][3] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][3];
            delamination_damage_affected_stress_matrix[i][4] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][4];
            delamination_damage_affected_stress_matrix[i][5] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][5];

            // if (negative_interfacial_stress_index [i] != 1.0 && negative_interfacial_stress_index [i+1] != 1.0) {
            //     delamination_damage_affected_stress_matrix[i][2] = (1-layer_damage_variable) * delamination_damage_affected_stress_matrix[i][2];
            // }
        }
        // Calculating output stresses

        // for(int i=0; i < mConstitutiveLaws.size(); ++i) {

        //     Properties& r_prop             = *(it_prop_begin + i);
        //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        //     const double factor            = mCombinationFactors[i];

        //     std::vector<double> layer_damage_vector(6);
        //     layer_damage_vector = {1,
        //                             1,
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1]),
        //                             1,
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1]),
        //                             (1-delamination_damage[i])*(1-delamination_damage[i+1])};

        //     auxiliar_stress_vector[0] += factor * layer_damage_vector[0] * layer_stress[i][0];
        //     auxiliar_stress_vector[1] += factor * layer_damage_vector[1] * layer_stress[i][1];
        //     auxiliar_stress_vector[2] += factor * layer_damage_vector[2] * layer_stress[i][2];
        //     auxiliar_stress_vector[3] += factor * layer_damage_vector[3] * layer_stress[i][3];
        //     auxiliar_stress_vector[4] += factor * layer_damage_vector[4] * layer_stress[i][4];
        //     auxiliar_stress_vector[5] += factor * layer_damage_vector[5] * layer_stress[i][5];
        // }

        // double damage_coeff = 1;

        for(int i=0; i < mConstitutiveLaws.size(); ++i) {
            const double factor = mCombinationFactors[i];
            delamination_damage_affected_stress_vector += factor * delamination_damage_affected_stress_matrix[i];
        }

        auxiliar_stress_vector = delamination_damage_affected_stress_vector;

        // auxiliar_stress_vector[0] = damage_coeff * undamaged_auxiliar_stress_vector[0];
        // auxiliar_stress_vector[1] = damage_coeff * undamaged_auxiliar_stress_vector[1];
        // auxiliar_stress_vector[2] = damage_coeff * undamaged_auxiliar_stress_vector[2];
        // auxiliar_stress_vector[3] = damage_coeff * undamaged_auxiliar_stress_vector[3];
        // auxiliar_stress_vector[4] = damage_coeff * undamaged_auxiliar_stress_vector[4];
        // auxiliar_stress_vector[5] = damage_coeff * undamaged_auxiliar_stress_vector[5];



        // End calculating output stresses



        // KRATOS_WATCH(status_coeff);
        
        noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

        

        // Previous flags restored
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");



    // End finalizing delamination damage

    //


    //

        // begin finaliaing delamination Damage Criterion V2

        // const double T0n = 16000000; // Interfacial Normal Strength
        // const double T0s = 27000000; // Interfacial Shear Strength
        // const double T0t = 27000000; // Interfacial Shear Strength
        // const double GIc = 102; // Mode I Energy Release Rate
        // const double GIIc = 194; // Mode II Energy Release Rate
        // const double Eta = 2; // Benzeggagh-Kenane (B-K) Law Coefficient
        // double Gc = mGc; // Mix Mode Energy Release Rate
        // double initial_threshold = minitial_threshold;
        // double threshold = mthreshold;
        // double delamination_damage = mdelamination_damage;
        // double AParameter = mAParameter;
        // double DamageIndicator = mDamageIndicator;
        // const double characteristic_length = 0.0001; // Characteristic Length of the Cohesive Part
        // const double tolerance = std::numeric_limits<double>::epsilon();
        // const double Fd = std::pow(undamaged_auxiliar_stress_vector[2]/T0n,2.0)+std::pow(undamaged_auxiliar_stress_vector[4]/T0s,2.0)+std::pow(undamaged_auxiliar_stress_vector[5]/T0t,2.0); // Damage Initiation Criterion
        // double T_eq = std::sqrt(std::pow(undamaged_auxiliar_stress_vector[2],2.0)+std::pow(undamaged_auxiliar_stress_vector[4],2.0)+std::pow(undamaged_auxiliar_stress_vector[5],2.0));
        // if (Fd >= 1.0 && DamageIndicator == 1) {
        //     initial_threshold = T_eq;
        //     mthreshold = T_eq;
        //     double Gn = undamaged_auxiliar_stress_vector[2] * strain_vector[2] / 2.0;
        //     double Gs = undamaged_auxiliar_stress_vector[4] * strain_vector[4] / 2.0;
        //     double Gt = undamaged_auxiliar_stress_vector[5] * strain_vector[5] / 2.0;
        //     double mode_mix_factor = (Gs+Gt) / (Gn+Gs+Gt);
        //     Gc = GIc + (GIIc - GIc) * std::pow(mode_mix_factor, Eta); // Benzeggagh-Kenane (B-K) Law
        //     // AParameter = -std::pow(initial_threshold, 2) / (2.0 * 800000000 * Gc / characteristic_length); // Linear
        //     AParameter = 1.00 / (Gc * 800000000 / (characteristic_length * std::pow(initial_threshold, 2)) - 0.5); // Exponential
        //     mGc = Gc;
        //     mAParameter = AParameter;
        //     minitial_threshold = initial_threshold;
        //     mDamageIndicator = 2;
        // }
        // const double Fp = T_eq - threshold;
        // if (Fp <= tolerance) { // Elastic case
        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5]; 
        // } else { // Damage case
        //     // delamination_damage = (1.0 - initial_threshold / T_eq) / (1.0 + AParameter); // Linear
        //     delamination_damage = 1.0 - (initial_threshold / T_eq) * std::exp(AParameter *
        //         (1.0 - T_eq / initial_threshold)); // Exponential

        //     auxiliar_stress_vector[0] = undamaged_auxiliar_stress_vector[0];
        //     auxiliar_stress_vector[1] = undamaged_auxiliar_stress_vector[1];
        //     auxiliar_stress_vector[2] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[2];
        //     auxiliar_stress_vector[3] = undamaged_auxiliar_stress_vector[3];
        //     auxiliar_stress_vector[4] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[4];
        //     auxiliar_stress_vector[5] = (1.0-delamination_damage) * undamaged_auxiliar_stress_vector[5];
        //     mthreshold = T_eq;
        //     mdelamination_damage = delamination_damage;
        // }

        // End finaliaing delamination Damage Criterion V2

    //


}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
    r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
    r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    // Previous flags saved
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
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
    r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
    r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::ResetMaterial(
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
void TractionSeparationLaw3D<TDim>::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int TractionSeparationLaw3D<TDim>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliar output
    int aux_out = 0;

    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "TractionSeparationLaw3D: No constitutive laws defined" << std::endl;

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
void TractionSeparationLaw3D<TDim>::CalculateRotationMatrix(
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

        BoundedMatrix<double, 3, 3>  rotation_matrix;

        if (std::abs(euler_angle_phi) + std::abs(euler_angle_theta) + std::abs(euler_angle_hi) > machine_tolerance) {
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(euler_angle_phi, euler_angle_theta, euler_angle_hi, rotation_matrix);
            AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(rotation_matrix, rRotationMatrix);
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
void TractionSeparationLaw3D<TDim>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues, const BaseType::StressMeasure& rStressMeasure)
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
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 4);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateAlmansiStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
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
void TractionSeparationLaw3D<TDim>::CalculateGreenLagrangeStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
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

template class TractionSeparationLaw3D<2>;
template class TractionSeparationLaw3D<3>;

} // Namespace Kratos
