// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Fernando Rastellini
//                   Tommaso Tudisco
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/composite_vanishing_fibre_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

CompositeVanishingFibreLaw::CompositeVanishingFibreLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

CompositeVanishingFibreLaw::CompositeVanishingFibreLaw(const std::vector<double>& rCombinationFactors) : ConstitutiveLaw()
{
    // We compute the proportion of the factors (must be over 1)
    double aux_factor = 0.0;
    for (IndexType i = 0; i < rCombinationFactors.size(); ++i) {
        aux_factor += rCombinationFactors[i];
    }

    KRATOS_ERROR_IF(aux_factor < std::numeric_limits<double>::epsilon()) << "Wrong factors in CompositeVanishingFibreLaw" << std::endl;

    // Resize
    mCombinationFactors.resize(rCombinationFactors.size());

    // We fill the maps
    for (IndexType i = 0; i < rCombinationFactors.size(); ++i) {
        mCombinationFactors[i] = rCombinationFactors[i]/aux_factor;
    }
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

CompositeVanishingFibreLaw::CompositeVanishingFibreLaw(const CompositeVanishingFibreLaw& rOther)
    : ConstitutiveLaw(rOther),
      mConstitutiveLaws(rOther.mConstitutiveLaws),  // TODO FRastellini: This is not ok.
      mCombinationFactors(rOther.mCombinationFactors)
{
    // TODO FRastellini: To perform a copy of the content of Components CLs
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer CompositeVanishingFibreLaw::Clone() const
{
    return Kratos::make_shared<CompositeVanishingFibreLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer CompositeVanishingFibreLaw::Create(Kratos::Parameters NewParameters) const
{
    // We do some checks
    KRATOS_ERROR_IF_NOT(NewParameters.Has("combination_factors")) << "CompositeVanishingFibreLaw: Please define combination_factors" << std::endl;

    const SizeType number_of_factors = NewParameters["combination_factors"].size();

    // We create the vectors
    std::vector<double> combination_factors(number_of_factors);

    for (IndexType i_layer = 0; i_layer < number_of_factors; ++i_layer) {
        combination_factors[i_layer] = NewParameters["combination_factors"][i_layer].GetDouble();
    }

    // We create the law
    return Kratos::make_shared<CompositeVanishingFibreLaw>(combination_factors);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

CompositeVanishingFibreLaw::~CompositeVanishingFibreLaw()
{
    // TODO FRastellini: To check. Is it necessary to destruct also the components CLs????
};

/***********************************************************************************/
/***********************************************************************************/

std::size_t CompositeVanishingFibreLaw::WorkingSpaceDimension()
{
    SizeType counter = 0;
    SizeType dimension = 0;
    // We perform the check in each material layer
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

std::size_t CompositeVanishingFibreLaw::GetStrainSize()
{
    SizeType counter = 0;
    SizeType strain_size = 0;
    // We perform the check in each material layer
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

bool CompositeVanishingFibreLaw::Has(const Variable<bool>& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<int>& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<double>& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<Vector>& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<Matrix>& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
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

bool CompositeVanishingFibreLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
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

bool& CompositeVanishingFibreLaw::GetValue(
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
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

int& CompositeVanishingFibreLaw::GetValue(
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
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

double& CompositeVanishingFibreLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    rValue = 0.0;
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        double aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

Vector& CompositeVanishingFibreLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        Vector aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& CompositeVanishingFibreLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        Matrix aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& CompositeVanishingFibreLaw::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(3);
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        array_1d<double, 3 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& CompositeVanishingFibreLaw::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(6);
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        array_1d<double, 6 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
    // WARNING FRastellini: Only the variable of the first found component is returned! Should we smooth the value of all components?
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::SetValue(
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

void CompositeVanishingFibreLaw::SetValue(
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

void CompositeVanishingFibreLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool& CompositeVanishingFibreLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer (for bools could be problematic)
    rValue = false;
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = material_properties.GetSubProperty(i + 1);

        rParameterValues.SetMaterialProperties(r_prop);
        bool aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        // rValue = rValue & aux_value; //and (intersection)
        rValue = rValue || aux_value; //or    (union)
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& CompositeVanishingFibreLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer (for integers could be problematic)
    rValue = 0;
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = material_properties.GetSubProperty(i + 1);

        rParameterValues.SetMaterialProperties(r_prop);
        int aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& CompositeVanishingFibreLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = 0.0;
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = material_properties.GetSubProperty(i + 1);

        rParameterValues.SetMaterialProperties(r_prop);
        double aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& CompositeVanishingFibreLaw::CalculateValue(
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
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
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

        rValue = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
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

        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else {
        const Properties& material_properties  = rParameterValues.GetMaterialProperties();

        // We combine the value of each layer
        rValue.clear();
        for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
            const double factor = mCombinationFactors[i];
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            Properties& r_prop = material_properties.GetSubProperty(i + 1);

            rParameterValues.SetMaterialProperties(r_prop);
            Vector aux_value;
            p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
            rValue += factor * aux_value;
        }

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& CompositeVanishingFibreLaw::CalculateValue(
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

        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false );
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

        rValue = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == DEFORMATION_GRADIENT) { // TODO: Make in the future modifications for take into account different layers combinations
        rValue = rParameterValues.GetDeformationGradientF();
    } else {
        const Properties& material_properties  = rParameterValues.GetMaterialProperties();

        // We combine the value of each layer
        rValue.clear();
        for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
            const double factor = mCombinationFactors[i];
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            Properties& r_prop = material_properties.GetSubProperty(i + 1);

            rParameterValues.SetMaterialProperties(r_prop);
            Matrix aux_value;
            p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
            rValue += factor * aux_value;
        }

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }

    return( rValue );
}
/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& CompositeVanishingFibreLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = ZeroVector(3);
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = material_properties.GetSubProperty(i + 1);

        rParameterValues.SetMaterialProperties(r_prop);
        array_1d<double, 3 > aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& CompositeVanishingFibreLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = ZeroVector(6);
    for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        const double factor = mCombinationFactors[i];
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = material_properties.GetSubProperty(i + 1);

        rParameterValues.SetMaterialProperties(r_prop);
        array_1d<double, 6 > aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

bool CompositeVanishingFibreLaw::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it component by component
    bool valid_input = true;
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        if (p_law->ValidateInput(r_prop)) {
            valid_input = false;
            break;
        }
    }

    return valid_input;
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StrainMeasure CompositeVanishingFibreLaw::GetStrainMeasure()
{
    // We return the first one
    for (auto& p_law : mConstitutiveLaws) {
        return p_law->GetStrainMeasure();
    }

    return StrainMeasure_Infinitesimal;
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure CompositeVanishingFibreLaw::GetStressMeasure()
{
    // We return the first one
    for (auto& p_law : mConstitutiveLaws) {
        return p_law->GetStressMeasure();
    }

    return StressMeasure_Cauchy;
}

/***********************************************************************************/
/***********************************************************************************/

bool CompositeVanishingFibreLaw::IsIncremental()
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

void CompositeVanishingFibreLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // Resizing first
    mConstitutiveLaws.resize(mCombinationFactors.size());

    // We create the inner constitutive laws
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);

        ConstitutiveLaw::Pointer p_inner_law = (r_prop)[CONSTITUTIVE_LAW]->Clone();
        p_inner_law->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
        mConstitutiveLaws[i] = p_inner_law;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the InitializeSolutionStep in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->InitializeSolutionStep(r_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the FinalizeSolutionStep in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->FinalizeSolutionStep(r_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the InitializeNonLinearIteration in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->InitializeNonLinearIteration(r_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the FinalizeNonLinearIteration in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->FinalizeNonLinearIteration(r_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
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

void  CompositeVanishingFibreLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same
    if(r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient;
        this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
        const Matrix C_matrix = prod(trans(F_deformation_gradient),F_deformation_gradient);
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

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (C_matrix - identity_matrix);
        // Green-Lagrangian Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
        Matrix constitutive_matrix = ZeroMatrix(voigt_size, voigt_size);
        for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
            Properties& r_prop = *(r_material_properties.GetSubProperties().begin() + i);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            const double factor = mCombinationFactors[i];

            rValues.SetMaterialProperties(r_prop);
            Matrix aux_value(voigt_size, voigt_size);
            p_law->CalculateValue(rValues, CONSTITUTIVE_MATRIX_PK2, aux_value);
            constitutive_matrix += factor * aux_value;
        }

        rValues.GetConstitutiveMatrix() = constitutive_matrix;
        rValues.SetMaterialProperties(r_material_properties);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector stress_vector = ZeroVector(voigt_size);
        for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
            Properties& r_prop = *(r_material_properties.GetSubProperties().begin() + i);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            const double factor = mCombinationFactors[i];

            rValues.SetMaterialProperties(r_prop);
            Vector aux_value(voigt_size);
            p_law->CalculateValue(rValues, PK2_STRESS_VECTOR, aux_value);
            stress_vector += factor * aux_value;
        }

        rValues.GetStressVector() = stress_vector;
        rValues.SetMaterialProperties(r_material_properties);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    if(r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient;
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
        Matrix inverse_B_tensor ( dimension, dimension );
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix( B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix constitutive_matrix = ZeroMatrix(voigt_size, voigt_size);
        for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
            Properties& r_prop = *(r_material_properties.GetSubProperties().begin() + i);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            const double factor = mCombinationFactors[i];

            rValues.SetMaterialProperties(r_prop);
            Matrix aux_value(voigt_size, voigt_size);
            p_law->CalculateValue(rValues, CONSTITUTIVE_MATRIX_KIRCHHOFF, aux_value);
            constitutive_matrix += factor * aux_value;
        }

        rValues.GetConstitutiveMatrix() = constitutive_matrix;
        rValues.SetMaterialProperties(r_material_properties);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector stress_vector = ZeroVector(voigt_size);
        for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
            Properties& r_prop = *(r_material_properties.GetSubProperties().begin() + i);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
            const double factor = mCombinationFactors[i];

            rValues.SetMaterialProperties(r_prop);
            Vector aux_value(voigt_size);
            p_law->CalculateValue(rValues, KIRCHHOFF_STRESS_VECTOR, aux_value);
            stress_vector += factor * aux_value;
        }

        rValues.GetStressVector() = stress_vector;
        rValues.SetMaterialProperties(r_material_properties);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
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

void CompositeVanishingFibreLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->InitializeMaterialResponsePK1(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->InitializeMaterialResponsePK2(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->InitializeMaterialResponseKirchhoff(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->InitializeMaterialResponseCauchy(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->FinalizeMaterialResponsePK1(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->FinalizeMaterialResponsePK2(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->FinalizeMaterialResponseKirchhoff(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = material_properties.GetSubProperty(i + 1);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        rValues.SetMaterialProperties(r_prop);
        p_law->FinalizeMaterialResponseCauchy(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We perform the reset in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        p_law->ResetMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    }
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/

void CompositeVanishingFibreLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

int CompositeVanishingFibreLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // The auxiliar output
    int aux_out = 0;

    // We perform the check in each layer
    for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];

        aux_out += p_law->Check(r_prop, rElementGeometry, rCurrentProcessInfo);
    }

    return aux_out;
}

} // Namespace Kratos
