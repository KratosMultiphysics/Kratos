//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_orthotropic_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticOrthotropic3DLaw::LinearElasticOrthotropic3DLaw()
    : HyperElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElasticOrthotropic3DLaw::LinearElasticOrthotropic3DLaw(const LinearElasticOrthotropic3DLaw& rOther)
    : HyperElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElasticOrthotropic3DLaw::Clone() const
{
    return Kratos::make_shared<LinearElasticOrthotropic3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElasticOrthotropic3DLaw::~LinearElasticOrthotropic3DLaw()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void  LinearElasticOrthotropic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//

    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
	//only needed
	const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();

        //4.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(DeformationGradientF),DeformationGradientF);

        //5.-Green-Lagrange Strain:

        //E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);
    }

    //7.-Calculate Total PK2 stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

      }
      else {

	Matrix ConstitutiveMatrix( StrainVector.size(), StrainVector.size() );
	noalias(ConstitutiveMatrix) = ZeroMatrix( StrainVector.size(), StrainVector.size() );

	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
      }

    }
    else if(  Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) && Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );

    }

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Strain "<<StrainVector<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;

}


//************************************************************************************
//************************************************************************************


void LinearElasticOrthotropic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();


    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//

    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(DeformationGradientF,trans(DeformationGradientF));

        //3.-Almansi Strain:

        // e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);

    }

    //7.-Calculate total Kirchhoff stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );

	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

      }
      else {

	Matrix ConstitutiveMatrix( StrainVector.size() ,StrainVector.size() );
	noalias(ConstitutiveMatrix) = ZeroMatrix( StrainVector.size() ,StrainVector.size() );

	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );

	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

      }

    }
    else if(  Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) && Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

      Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
      this->CalculateLinearElasticMatrix( ConstitutiveMatrix, MaterialProperties );

    }

    //std::cout<<" Strain "<<StrainVector<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;

}



//***********************COMPUTE TOTAL STRESS PK2*************************************
//************************************************************************************


void LinearElasticOrthotropic3DLaw::CalculateStress( const Vector & rStrainVector,
        const Matrix & rConstitutiveMatrix,
        Vector& rStressVector )
{

    //1.-2nd Piola Kirchhoff StressVector increment
    if( rStressVector.size() != rStrainVector.size() )
      rStressVector.resize(rStrainVector.size(),false);

    noalias(rStressVector) = prod(rConstitutiveMatrix,rStrainVector);


}



//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************


void LinearElasticOrthotropic3DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
																  const Properties& rMaterialProperties )
{
	double E1 = rMaterialProperties[YOUNG_MODULUS_X];
	double E2 = rMaterialProperties[YOUNG_MODULUS_Y];
	double E3 = rMaterialProperties[YOUNG_MODULUS_Z];

	double v12 = rMaterialProperties[POISSON_RATIO_XY];
	double v23 = rMaterialProperties[POISSON_RATIO_YZ];
	double v13 = rMaterialProperties[POISSON_RATIO_XZ];

	double P1 = 1.0/(E2*E2*v12*v12 + 2.0*E3*E2*v12*v13*v23 + E3*E2*v13*v13 - E1*E2 + E1*E3*v23*v23);
	double P2 = E1*E1;
	double P3 = E2*E2;
	double P4 = E1*v23 + E2*v12*v13;
	double P5 = E2*v12 + E3*v13*v23;
	double P6 = E3*E3;

    rConstitutiveMatrix.clear();

	rConstitutiveMatrix(0, 0) = -P1*P2*(- E3*v23*v23 + E2);
	rConstitutiveMatrix(0, 1) = -E1*E2*P1*P5;
	rConstitutiveMatrix(0, 2) = -E2*E3*P1*(E1*v13 + E1*v12*v23);
	rConstitutiveMatrix(1, 0) = -E1*E2*P1*P5;
	rConstitutiveMatrix(1, 1) = -P1*P3*(- E3*v13*v13 + E1);
	rConstitutiveMatrix(1, 2) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 0) = -E1*E2*E3*P1*(v13 + v12*v23);
	rConstitutiveMatrix(2, 1) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 2) = -E2*E3*P1*(- E2*v12*v12 + E1);
	rConstitutiveMatrix(3, 3) = (E2*P2)/(P2 + v12*(P2 + P3) + E1*E2)/2.0;
	rConstitutiveMatrix(4, 4) = (E3*P3)/(P3 + v23*(P3 + P6) + E2*E3)/2.0;
	rConstitutiveMatrix(5, 5) = (E3*P2)/(P2 + v13*(P2 + P6) + E1*E3)/2.0;
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElasticOrthotropic3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ANISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

bool LinearElasticOrthotropic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int LinearElasticOrthotropic3DLaw::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS_X.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_X))
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS_X has Key zero or invalid value ", "" )

	if(YOUNG_MODULUS_Y.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Y))
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS_Y has Key zero or invalid value ", "" )

	if(YOUNG_MODULUS_Z.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Z))
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS_Z has Key zero or invalid value ", "" )

    if(POISSON_RATIO_XY.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XY))
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO_XY has Key zero invalid value ", "" )

	if(POISSON_RATIO_YZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_YZ))
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO_YZ has Key zero invalid value ", "" )

	if(POISSON_RATIO_XZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XZ))
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO_XZ has Key zero invalid value ", "" )

    if(DENSITY.Key() == 0 || !rMaterialProperties.Has(DENSITY))
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero or invalid value ", "" )


    return 0;

}


} // Namespace Kratos
