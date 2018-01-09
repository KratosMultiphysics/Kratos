// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes
// #include<cmath>

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropic3D::ElasticIsotropic3D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ElasticIsotropic3D::ElasticIsotropic3D(const ElasticIsotropic3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticIsotropic3D::Clone() const
{
    ElasticIsotropic3D::Pointer p_clone(new ElasticIsotropic3D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropic3D::~ElasticIsotropic3D()
{
};

//************************************************************************************
//************************************************************************************

void  ElasticIsotropic3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    const double& E          = MaterialProperties[YOUNG_MODULUS];
    const double& NU    = MaterialProperties[POISSON_RATIO];

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        CalculateCauchyGreenStrain(rValues, StrainVector);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( ConstitutiveMatrix, E, NU );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        if (rValues.IsSetDeformationGradientF() == true)
        {
            CalculateCauchyGreenStrain(rValues, StrainVector);
        }

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {
            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            noalias(StressVector) = prod(ConstitutiveMatrix, StrainVector);
        }
        else
        {
            CalculatePK2Stress( StrainVector, StressVector, E, NU );
        }
    }
}

//************************************************************************************
//************************************************************************************

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything

void ElasticIsotropic3D::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    // TODO: Add if necessary
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // TODO: Add if necessary
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // TODO: Add if necessary
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    // TODO: Add if necessary
}

//************************************************************************************
//************************************************************************************

double& ElasticIsotropic3D::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    const Properties& MaterialProperties  = rParameterValues.GetMaterialProperties();
    Vector& StrainVector                  = rParameterValues.GetStrainVector();
    Vector& StressVector                  = rParameterValues.GetStressVector();
    const double& E          = MaterialProperties[YOUNG_MODULUS];
    const double& NU    = MaterialProperties[POISSON_RATIO];
    
    if (rThisVariable == STRAIN_ENERGY)
    {
        CalculateCauchyGreenStrain(rParameterValues, StrainVector);
        CalculatePK2Stress( StrainVector, StressVector, E, NU );

        rValue = 0.5 * inner_prod(StrainVector,StressVector); // Strain energy = 0.5*E:C:E
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void ElasticIsotropic3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 6;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int ElasticIsotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.0)
    {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
    {
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
    }

    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.0)
    {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }

    return 0;
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::CalculateElasticMatrix(
    Matrix& C,
    const double E,
    const double NU
)
{
    const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    const double c2 = c1 * ( 1 - NU );
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 0, 3 ) = 0.0;
    C( 0, 4 ) = 0.0;
    C( 0, 5 ) = 0.0;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 1, 3 ) = 0.0;
    C( 1, 4 ) = 0.0;
    C( 1, 5 ) = 0.0;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 2, 3 ) = 0.0;
    C( 2, 4 ) = 0.0;
    C( 2, 5 ) = 0.0;
    C( 3, 0 ) = 0.0;
    C( 3, 1 ) = 0.0;
    C( 3, 2 ) = 0.0;
    C( 3, 3 ) = c4;
    C( 3, 4 ) = 0.0;
    C( 3, 5 ) = 0.0;
    C( 4, 0 ) = 0.0;
    C( 4, 1 ) = 0.0;
    C( 4, 2 ) = 0.0;
    C( 4, 3 ) = 0.0;
    C( 4, 4 ) = c4;
    C( 4, 5 ) = 0.0;
    C( 5, 0 ) = 0.0;
    C( 5, 1 ) = 0.0;
    C( 5, 2 ) = 0.0;
    C( 5, 3 ) = 0.0;
    C( 5, 4 ) = 0.0;
    C( 5, 5 ) = c4;
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    const double E,
    const double NU
)
{
    SizeType SizeSystem = GetStrainSize();
    Matrix C(SizeSystem,SizeSystem);
    CalculateElasticMatrix(C, E, NU);
    noalias(rStressVector) = prod(C,rStrainVector);
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropic3D::CalculateCauchyGreenStrain(
    Parameters& rValues,
    Vector& rStrainVector
)
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    Matrix Etensor = prod(trans(F),F);
    Etensor -= IdentityMatrix(3,3);
    Etensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(Etensor);
}

} // Namespace Kratos
