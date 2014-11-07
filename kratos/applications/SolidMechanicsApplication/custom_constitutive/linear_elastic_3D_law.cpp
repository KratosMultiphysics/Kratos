//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/linear_elastic_3D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElastic3DLaw::LinearElastic3DLaw()
    : HyperElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElastic3DLaw::LinearElastic3DLaw(const LinearElastic3DLaw& rOther)
    : HyperElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElastic3DLaw::Clone() const
{
    LinearElastic3DLaw::Pointer p_clone(new LinearElastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElastic3DLaw::~LinearElastic3DLaw()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void  LinearElastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
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

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {

	//only needed 
	const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
		
	//2.-Total Deformation Gradient
        Matrix TotalDeformationGradientF0 = DeformationGradientF;


        //4.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(TotalDeformationGradientF0),TotalDeformationGradientF0);

        //5.-Green-Lagrange Strain:

        //E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);


	if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION )  || Options.Is( ConstitutiveLaw::FINAL_CONFIGURATION )){

	  //2.-Total Deformation Gradient
	  Matrix& DeformationGradientF0   = rValues.GetDeformationGradientF0();

	  TotalDeformationGradientF0      = prod(DeformationGradientF, DeformationGradientF0);
	  
	  //3.-Left Cauchy Green
	  Matrix LeftCauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));
	  
	  //4.-Almansi Strain
	  // e= 0.5*(1-invbT*invb)
	  this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);
 
	  //5.-Pull-back Almansi Strain
	  this->TransformStrains(StrainVector,DeformationGradientF,StrainMeasure_Almansi,StrainMeasure_GreenLagrange);
	  
	}


      }

    //7.-Calculate Total PK2 stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
      	
	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

      }
      else {

	Matrix ConstitutiveMatrix = ZeroMatrix( StrainVector.size() );
	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
      }
      
    }
    else if(  Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) && Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

    }

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Strain "<<StrainVector<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;

}


//************************************************************************************
//************************************************************************************


void LinearElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
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

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	//1.-Compute total deformation gradient
        Matrix& DeformationGradientF0      = rValues.GetDeformationGradientF0();

	Matrix TotalDeformationGradientF0  = prod(DeformationGradientF, DeformationGradientF0);
	TotalDeformationGradientF0         = this->DeformationGradient3D( TotalDeformationGradientF0 );

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));

        //3.-Almansi Strain:

        // e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);

      }

    //7.-Calculate total Kirchhoff stress

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
      
      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
	
	Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
      
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

      }
      else {
	
	Matrix ConstitutiveMatrix = ZeroMatrix( StrainVector.size() );
	
	this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
      
	this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
      
      }
      
    }
    else if(  Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) && Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

      Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
      this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

    }

    //std::cout<<" Strain "<<StrainVector<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;

}



//***********************COMPUTE TOTAL STRESS PK2*************************************
//************************************************************************************


void LinearElastic3DLaw::CalculateStress( const Vector & rStrainVector,
        const Matrix & rConstitutiveMatrix,
        Vector& rStressVector )
{

    //1.-2nd Piola Kirchhoff StressVector increment
    rStressVector = prod(rConstitutiveMatrix,rStrainVector);


}



//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************


void LinearElastic3DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
        const double &rYoungModulus,
        const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    //plane strain constitutive matrix: 
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1-2*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
    rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
    rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );

}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElastic3DLaw::GetLawFeatures(Features& rFeatures)
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

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

bool LinearElastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int LinearElastic3DLaw::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
        KRATOS_ERROR( std::invalid_argument,"DENSITY has Key zero or invalid value ", "" )


    return 0;

}


} // Namespace Kratos
