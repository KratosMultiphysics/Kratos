// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_constitutive/user_provided_linear_elastic_law.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
UserProvidedLinearElasticLaw<TDim>::UserProvidedLinearElasticLaw()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template<unsigned int TDim>
UserProvidedLinearElasticLaw<TDim>::UserProvidedLinearElasticLaw(const UserProvidedLinearElasticLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::Pointer UserProvidedLinearElasticLaw<TDim>::Clone() const
{
    return Kratos::make_shared<UserProvidedLinearElasticLaw<TDim>>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
UserProvidedLinearElasticLaw<TDim>::~UserProvidedLinearElasticLaw()
{
};

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  UserProvidedLinearElasticLaw<TDim>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get the constitutive law options
    Flags & r_constitutive_law_options = rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if(r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrainVector(rValues, r_strain_vector);
    }

    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
    }

    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("");
}

// /***********************************************************************************/
// /***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything
template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

// /***********************************************************************************/
// /***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything
template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

// /***********************************************************************************/
// /***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything
template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& UserProvidedLinearElasticLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    Vector& r_strain_vector = rParameterValues.GetStrainVector();
    Vector& r_stress_vector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        this->CalculateGreenLagrangeStrainVector(rParameterValues, r_strain_vector);
        this->CalculatePK2Stress( r_strain_vector, r_stress_vector, rParameterValues);
        rValue = 0.5 * inner_prod( r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& UserProvidedLinearElasticLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN || rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR || rThisVariable == ALMANSI_STRAIN_VECTOR) {
        this->CalculateGreenLagrangeStrainVector( rParameterValues, rValue);
    } else if (rThisVariable == STRESSES || rThisVariable == CAUCHY_STRESS_VECTOR || rThisVariable == KIRCHHOFF_STRESS_VECTOR || rThisVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // We compute the stress
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );
        UserProvidedLinearElasticLaw<TDim>::CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Matrix& UserProvidedLinearElasticLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX || rThisVariable == CONSTITUTIVE_MATRIX_PK2 || rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int UserProvidedLinearElasticLaw<TDim>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK_VARIABLE_KEY(ELASTICITY_TENSOR);
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    noalias(rConstitutiveMatrix) = r_material_properties[ELASTICITY_TENSOR];
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const Matrix C = r_material_properties[ELASTICITY_TENSOR];
    noalias(rStressVector) = prod(C, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void UserProvidedLinearElasticLaw<TDim>::CalculateGreenLagrangeStrainVector(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    const SizeType dim = this->WorkingSpaceDimension();

    // Get the deformation gradient tensor F
    const Matrix& rF = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(rF.size1()!= dim || rF.size2() != dim) << "expected size of F " << dim << "x" << dim << ", got " << rF.size1() << "x" << rF.size2() << std::endl;

    // Calculate the Cauchy - Green strain tensor
    Matrix left_cauchy_green = prod(trans(rF), rF);

    // Calculate Green - Lagrange strain tensor
    ConstitutiveLawUtilities<StrainSize>::CalculateGreenLagrangianStrain(left_cauchy_green, rStrainVector);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class UserProvidedLinearElasticLaw<2>;
template class UserProvidedLinearElasticLaw<3>;

} // Namespace Kratos
