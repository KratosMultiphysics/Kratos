// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

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
    // We fill the maps
    for (IndexType i = 0; i < rSubPropertiesIDs.size(); ++i) {
        const IndexType id = rSubPropertiesIDs[i];
        mCombinationFactors.insert(std::pair<IndexType, double>({id, rCombinationFactors[i]}));
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

void RuleOfMixturesLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    for (auto& factors : mCombinationFactors) {
//         const Properties::Pointer p_prop = material_properties.GetSubProperty(factors.first);


    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f          = rValues.GetDeterminantF();

    TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

void  RuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    
    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const SizeType dimension = WorkingSpaceDimension();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the right Cauchy-Green tensor (C):
    const Matrix C_tensor = prod(trans( deformation_gradient_f), deformation_gradient_f);

    // Inverse of the right Cauchy-Green tensor (C):
    double aux_det;
    Matrix inverse_C_tensor(dimension, dimension);
    MathUtils<double>::InvertMatrix( C_tensor, inverse_C_tensor, aux_det);

    if(r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
//         CalculateCauchyGreenStrain(rValues, strain_vector);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
//         CalculateConstitutiveMatrixPK2( constitutive_matrix, inverse_C_tensor, determinant_f, lame_lambda, lame_mu );
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector& stress_vector = rValues.GetStressVector();
//         CalculatePK2Stress( inverse_C_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags& r_flags=rValues.GetOptions();

    const Properties& material_properties  = rValues.GetMaterialProperties();
    Vector& strain_vector                  = rValues.GetStrainVector();
    Vector& stress_vector                  = rValues.GetStressVector();

    // The material properties
    const double young_modulus = material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = material_properties[POISSON_RATIO];

    // The deformation gradient
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f = rValues.GetDeterminantF();
    KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;

    // The LAME parameters
    const double lame_lambda = (young_modulus * poisson_coefficient)/((1.0 + poisson_coefficient)*(1.0 - 2.0 * poisson_coefficient));
    const double lame_mu = young_modulus/(2.0 * (1.0 + poisson_coefficient));

    // We compute the left Cauchy-Green tensor (B):
    const Matrix B_tensor = prod(deformation_gradient_f, trans( deformation_gradient_f));

    if(r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
//         CalculateAlmansiStrain(rValues, strain_vector);
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
//         CalculateConstitutiveMatrixKirchhoff( constitutive_matrix, determinant_f, lame_lambda, lame_mu );
    }

    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
//         CalculateKirchhoffStress( B_tensor, stress_vector, determinant_f, lame_lambda, lame_mu );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseKirchhoff(rValues);

    Vector& stress_vector       = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    stress_vector       /= determinant_f;
    constitutive_matrix /= determinant_f;
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK1(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponsePK2(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseCauchy(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
}

/***********************************************************************************/
/***********************************************************************************/

void RuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
//     rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
//     this->CalculateMaterialResponseKirchhoff(rValues);
//     rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
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
    double value = 0.0;
    for (auto& factors : mCombinationFactors) {
//         Properties::Pointer p_prop = material_properties.GetSubProperty(factors.first);
        const double factor = factors.second;
    }

    // Reset properties
    rParameterValues.SetMaterialProperties(material_properties);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& RuleOfMixturesLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    ) 
{
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
        
        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
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
        
        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
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
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    const double nu = rMaterialProperties[POISSON_RATIO];
    const bool check = static_cast<bool>((nu >0.499 && nu<0.501) || (nu < -0.999 && nu > -1.01));
    KRATOS_ERROR_IF(check) << "POISSON_RATIO is invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is invalid value " << std::endl;

    return 0;
}

} // Namespace Kratos
