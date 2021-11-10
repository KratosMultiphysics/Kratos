// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_K0_3d_law.h"
#include "includes/checks.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ElasticIsotropicK03DLaw::ElasticIsotropicK03DLaw()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ElasticIsotropicK03DLaw::ElasticIsotropicK03DLaw(const ElasticIsotropicK03DLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ElasticIsotropicK03DLaw::Clone() const
{
    return Kratos::make_shared<ElasticIsotropicK03DLaw>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

ElasticIsotropicK03DLaw::~ElasticIsotropicK03DLaw()
{
}

/***********************************************************************************/
/***********************************************************************************/

void  ElasticIsotropicK03DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    // b.- Get Values to compute the constitutive law:
    Flags & r_options=rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }

    if( r_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( r_constitutive_matrix, rValues);
    }

    if( r_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
    }
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything

void ElasticIsotropicK03DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

double& ElasticIsotropicK03DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                const Variable<double>& rThisVariable,
                                                double& rValue)
{
    Vector& r_strain_vector = rParameterValues.GetStrainVector();
    Vector& r_stress_vector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        this->CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        this->CalculatePK2Stress( r_strain_vector, r_stress_vector, rParameterValues);

        rValue = 0.5 * inner_prod( r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ElasticIsotropicK03DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                const Variable<Vector>& rThisVariable,
                                                Vector& rValue)
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
        this->CalculateCauchyGreenStrain( rParameterValues, rValue);
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
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // We compute the stress
        ElasticIsotropicK03DLaw::CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain );
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

    return( rValue );
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ElasticIsotropicK03DLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                const Variable<Matrix>& rThisVariable,
                                                Matrix& rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void ElasticIsotropicK03DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}


void ElasticIsotropicK03DLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {}


/***********************************************************************************/
/***********************************************************************************/

int ElasticIsotropicK03DLaw::Check(const Properties& rMaterialProperties,
                                   const GeometryType& rElementGeometry,
                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(!rMaterialProperties.Has(YOUNG_MODULUS))
                    << "YOUNG_MODULUS is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is invalid value " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(POISSON_RATIO))
                    << "POISSON_RATIO is not availabe in material parameters" << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = ((nu >0.499 && nu<0.501) || (nu < -0.999 && nu > -1.01));
    KRATOS_ERROR_IF(check) << "POISSON_RATIO has invalid value " << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CheckClearElasticMatrix(Matrix& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CalculateElasticMatrix(Matrix& rConstitutiveMatrix,
                                                     ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(rConstitutiveMatrix);

    const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    const double c2 = c1 * ( 1 - NU );
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

    rConstitutiveMatrix( 0, 0 ) = c2;
    rConstitutiveMatrix( 0, 1 ) = c3;
    rConstitutiveMatrix( 0, 2 ) = c3;
    rConstitutiveMatrix( 1, 0 ) = c3;
    rConstitutiveMatrix( 1, 1 ) = c2;
    rConstitutiveMatrix( 1, 2 ) = c3;
    rConstitutiveMatrix( 2, 0 ) = c3;
    rConstitutiveMatrix( 2, 1 ) = c3;
    rConstitutiveMatrix( 2, 2 ) = c2;
    rConstitutiveMatrix( 3, 3 ) = c4;
    rConstitutiveMatrix( 4, 4 ) = c4;
    rConstitutiveMatrix( 5, 5 ) = c4;
}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                                 Vector& rStressVector,
                                                 ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c1 = E / ((1.00 + NU) * (1 - 2 * NU));
    const double c2 = c1 * (1 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1 - 2 * NU);

    rStressVector[0] = c2 * rStrainVector[0] + c3 * rStrainVector[1] + c3 * rStrainVector[2];
    rStressVector[1] = c3 * rStrainVector[0] + c2 * rStrainVector[1] + c3 * rStrainVector[2];
    rStressVector[2] = c3 * rStrainVector[0] + c3 * rStrainVector[1] + c2 * rStrainVector[2];
    rStressVector[3] = c4 * rStrainVector[3];
    rStressVector[4] = c4 * rStrainVector[4];
    rStressVector[5] = c4 * rStrainVector[5];

    // apply K0 procedure:
    const double& K0ValueXX = r_material_properties[K0_VALUE_XX];
    const double& K0ValueYY = r_material_properties[K0_VALUE_YY];
    const double& K0ValueZZ = r_material_properties[K0_VALUE_ZZ];
    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];

    if (K0MainDirection == VOIGT_INDEX_XX)
    {
        rStressVector[VOIGT_INDEX_YY] = K0ValueYY * rStressVector[VOIGT_INDEX_XX];
        rStressVector[VOIGT_INDEX_ZZ] = K0ValueZZ * rStressVector[VOIGT_INDEX_XX];
    }
    else if (K0MainDirection == VOIGT_INDEX_YY)
    {
        rStressVector[VOIGT_INDEX_XX] = K0ValueXX * rStressVector[VOIGT_INDEX_YY];
        rStressVector[VOIGT_INDEX_ZZ] = K0ValueZZ * rStressVector[VOIGT_INDEX_YY];
    }
    else if (K0MainDirection == VOIGT_INDEX_ZZ)
    {
        rStressVector[VOIGT_INDEX_XX] = K0ValueXX * rStressVector[VOIGT_INDEX_ZZ];
        rStressVector[VOIGT_INDEX_YY] = K0ValueYY * rStressVector[VOIGT_INDEX_ZZ];
    }
    else
    {
         KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticK03DLaw: " << K0MainDirection << std::endl;
    }

}

/***********************************************************************************/
/***********************************************************************************/

void ElasticIsotropicK03DLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector)
{
    const SizeType space_dimension = this->WorkingSpaceDimension();

    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();
    KRATOS_DEBUG_ERROR_IF(F.size1()!= space_dimension || F.size2() != space_dimension)
        << "expected size of F " << space_dimension << "x" << space_dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

    Matrix E_tensor = prod(trans(F),F);
    for(unsigned int i=0; i<space_dimension; ++i)
      E_tensor(i,i) -= 1.0;
    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos
