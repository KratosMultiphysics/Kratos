// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Fernando Rastellini
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/rule_of_mixtures_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

RuleOfMixturesLaw::RuleOfMixturesLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

RuleOfMixturesLaw::RuleOfMixturesLaw(
    const std::vector<IndexType>& rSubPropertiesIDs,
    const std::vector<double>& rCombinationFactors,
    const std::vector<double>& rMaterialRotationAngles
    ) : ConstitutiveLaw()
{
    // We compute the proportion of the factors (must be over 1)
    double aux_factor = 0.0;
    for (IndexType i = 0; i < rSubPropertiesIDs.size(); ++i) {
        aux_factor += rCombinationFactors[i];
    }

    KRATOS_ERROR_IF(aux_factor < std::numeric_limits<double>::epsilon()) << "Wrong factors in RuleOfMixturesLaw" << std::endl;

    // We fill the maps
    for (IndexType i = 0; i < rSubPropertiesIDs.size(); ++i) {
        const IndexType id = rSubPropertiesIDs[i];
        mCombinationFactors.insert(std::pair<IndexType, double>({id, rCombinationFactors[i]/aux_factor}));
        mMaterialRotationAngles.insert(std::pair<IndexType, double>({id, rMaterialRotationAngles[i]}));
    }
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

RuleOfMixturesLaw::RuleOfMixturesLaw(const RuleOfMixturesLaw& rOther)
    : ConstitutiveLaw(rOther),
      mCombinationFactors(rOther.mCombinationFactors),
      mMaterialRotationAngles(rOther.mMaterialRotationAngles)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer RuleOfMixturesLaw::Clone() const
{
    return Kratos::make_shared<RuleOfMixturesLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer RuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{
    // We do some checks
    KRATOS_ERROR_IF_NOT(NewParameters.Has("sub_properties_indexes")) << "RuleOfMixturesLaw: Please define sub_properties_indexes" << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("combination_factors")) << "RuleOfMixturesLaw: Please define combination_factors" << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("material_rotation_angles")) << "RuleOfMixturesLaw: Please define material_rotation_angles" << std::endl;

    const SizeType number_of_layers = NewParameters["sub_properties_indexes"].size();
    const SizeType number_of_factors = NewParameters["combination_factors"].size();
    const SizeType number_of_angles = NewParameters["material_rotation_angles"].size();

    KRATOS_ERROR_IF(number_of_layers != number_of_factors) << "The vectors sub_properties_indexes and combination_factors must have the same size" << std::endl;
    KRATOS_ERROR_IF(number_of_layers != number_of_angles) << "The vectors sub_properties_indexes and material_rotation_angles must have the same size" << std::endl;

    // We create the vectors
    std::vector<IndexType> sub_properties_ids(number_of_layers);
    std::vector<double> combination_factors(number_of_layers);
    std::vector<double> rotation_angles(number_of_layers);

    for (IndexType i_layer = 0; i_layer < number_of_layers; ++i_layer) {
        sub_properties_ids[i_layer] = NewParameters["sub_properties_indexes"][i_layer].GetInt();
        combination_factors[i_layer] = NewParameters["combination_factors"][i_layer].GetDouble();
        rotation_angles[i_layer] = NewParameters["material_rotation_angles"][i_layer].GetDouble();
    }

    // We create the law
    return Kratos::make_shared<RuleOfMixturesLaw>(sub_properties_ids, combination_factors, rotation_angles);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

RuleOfMixturesLaw::~RuleOfMixturesLaw()
{
};

/***********************************************************************************/
/***********************************************************************************/

std::size_t RuleOfMixturesLaw::WorkingSpaceDimension()
{
    SizeType counter = 0;
    SizeType dimension = 0;
    // We perform the check in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (counter == 0) {
            dimension = p_law->WorkingSpaceDimension();
        } else {
            KRATOS_ERROR_IF_NOT(dimension == p_law->WorkingSpaceDimension()) << "Combinig different size laws" << std::endl;
        }

        ++counter;
    }

    return dimension;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t RuleOfMixturesLaw::GetStrainSize()
{
    SizeType counter = 0;
    SizeType strain_size = 0;
    // We perform the check in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (counter == 0) {
            strain_size = p_law->GetStrainSize();
        } else {
            KRATOS_ERROR_IF_NOT(strain_size == p_law->GetStrainSize()) << "Combinig different size laws" << std::endl;
        }

        ++counter;
    }

    return strain_size;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<bool>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<int>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<double>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<Vector>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<Matrix>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            has = true;
            break;
        }
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool& RuleOfMixturesLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    // At least one layer should have the value
    rValue = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->GetValue(rThisVariable, rValue))
            break;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& RuleOfMixturesLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    // At least one layer should have the value
    rValue = 0;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        if (p_law->Has(rThisVariable)) {
            p_law->GetValue(rThisVariable, rValue);
            break;
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& RuleOfMixturesLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    rValue = 0.0;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        double aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& RuleOfMixturesLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        Vector aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& RuleOfMixturesLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        Matrix aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& RuleOfMixturesLaw::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(3);
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        array_1d<double, 3 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& RuleOfMixturesLaw::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(6);
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        array_1d<double, 6 > aux_value;
        p_law->GetValue(rThisVariable, aux_value);
        rValue += aux_value * factor;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        const double factor = factors.second;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->SetValue(rThisVariable, factor * rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool& RuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer (for bools could be problematic)
    rValue = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        const double factor = factors.second;

        rParameterValues.SetMaterialProperties(*p_prop);
        bool aux_value;
        p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        rValue += factor * aux_value;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& RuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer (for integers could be problematic)
    rValue = 0;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        const double factor = factors.second;

        rParameterValues.SetMaterialProperties(*p_prop);
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

double& RuleOfMixturesLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = 0.0;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        const double factor = factors.second;

        rParameterValues.SetMaterialProperties(*p_prop);
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

Vector& RuleOfMixturesLaw::CalculateValue(
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
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rParameterValues.SetMaterialProperties(*p_prop);
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

Matrix& RuleOfMixturesLaw::CalculateValue(
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
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rParameterValues.SetMaterialProperties(*p_prop);
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

array_1d<double, 3 >& RuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = ZeroVector(3);
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        const double factor = factors.second;

        rParameterValues.SetMaterialProperties(*p_prop);
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

array_1d<double, 6 >& RuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    const Properties& material_properties  = rParameterValues.GetMaterialProperties();

    // We combine the value of each layer
    rValue = ZeroVector(6);
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        const double factor = factors.second;

        rParameterValues.SetMaterialProperties(*p_prop);
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

bool RuleOfMixturesLaw::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it layer by layer
    bool valid_input = true;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        if (p_law->ValidateInput(*p_prop)) {
            valid_input = false;
            break;
        }
    }

    return valid_input;
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StrainMeasure RuleOfMixturesLaw::GetStrainMeasure()
{
    // We return the first one
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        return p_law->GetStrainMeasure();
    }
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure RuleOfMixturesLaw::GetStressMeasure()
{
    // We return the first one
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        return p_law->GetStressMeasure();
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool RuleOfMixturesLaw::IsIncremental()
{
    // We check it layer by layer
    bool is_incremental = false;
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
        if (p_law->IsIncremental()) {
            is_incremental = true;
            break;
        }
    }

    return is_incremental;
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We create the inner constitutive laws
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);

        ConstitutiveLaw::Pointer p_inner_law = (*p_prop)[CONSTITUTIVE_LAW]->Clone();
        p_inner_law->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
        mConstitutiveLaws.insert(std::pair<IndexType, ConstitutiveLaw::Pointer>({id, p_inner_law}));
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the InitializeSolutionStep in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->InitializeSolutionStep(*p_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the FinalizeSolutionStep in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->FinalizeSolutionStep(*p_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the InitializeNonLinearIteration in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->InitializeNonLinearIteration(*p_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We perform the FinalizeNonLinearIteration in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->FinalizeNonLinearIteration(*p_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
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

void  RuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
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
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = r_material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rValues.SetMaterialProperties(*p_prop);
            Matrix aux_value(voigt_size, voigt_size);
            p_law->CalculateValue(rValues, CONSTITUTIVE_MATRIX_PK2, aux_value);
            constitutive_matrix += factor * aux_value;
        }

        rValues.GetConstitutiveMatrix() = constitutive_matrix;
        rValues.SetMaterialProperties(r_material_properties);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector stress_vector = ZeroVector(voigt_size);
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = r_material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rValues.SetMaterialProperties(*p_prop);
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

void RuleOfMixturesLaw::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
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
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = r_material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rValues.SetMaterialProperties(*p_prop);
            Matrix aux_value(voigt_size, voigt_size);
            p_law->CalculateValue(rValues, CONSTITUTIVE_MATRIX_KIRCHHOFF, aux_value);
            constitutive_matrix += factor * aux_value;
        }

        rValues.GetConstitutiveMatrix() = constitutive_matrix;
        rValues.SetMaterialProperties(r_material_properties);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector stress_vector = ZeroVector(voigt_size);
        for (auto& factors : mCombinationFactors) {
            const IndexType id = factors.first;
            Properties::Pointer p_prop = r_material_properties.GetSubProperty(id);
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];
            const double factor = factors.second;

            rValues.SetMaterialProperties(*p_prop);
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

void RuleOfMixturesLaw::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
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

void RuleOfMixturesLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->InitializeMaterialResponsePK1(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->InitializeMaterialResponsePK2(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->InitializeMaterialResponseKirchhoff(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->InitializeMaterialResponseCauchy(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->FinalizeMaterialResponsePK1(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->FinalizeMaterialResponsePK2(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->FinalizeMaterialResponseKirchhoff(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& material_properties = rValues.GetMaterialProperties();

    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = material_properties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        rValues.SetMaterialProperties(*p_prop);
        p_law->FinalizeMaterialResponseCauchy(rValues);
    }

    rValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We perform the reset in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        p_law->ResetMaterial(*p_prop, rElementGeometry, rShapeFunctionsValues);
    }
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/

void RuleOfMixturesLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

int RuleOfMixturesLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // The auxiliar output
    int aux_out = 0;

    // We perform the check in each layer
    for (auto& factors : mCombinationFactors) {
        const IndexType id = factors.first;
        Properties::Pointer p_prop = rMaterialProperties.GetSubProperty(id);
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[id];

        aux_out += p_law->Check(*p_prop, rElementGeometry, rCurrentProcessInfo);
    }

    return aux_out;
}

} // Namespace Kratos
