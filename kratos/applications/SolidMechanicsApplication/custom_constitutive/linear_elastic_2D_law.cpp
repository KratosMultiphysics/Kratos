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
#include "custom_constitutive/linear_elastic_2D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LinearElastic2DLaw::LinearElastic2DLaw()
  : ConstitutiveLaw()
  {
    //pointer to constitutive law in properties :: serializer has nan in not set internal variables
    mDetF0 = 0;

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LinearElastic2DLaw::LinearElastic2DLaw(const LinearElastic2DLaw& rOther)
  : ConstitutiveLaw()
    ,mStressVector(rOther.mStressVector)
    ,mStressMeasure(rOther.mStressMeasure)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDetF0(rOther.mDetF0)
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer LinearElastic2DLaw::Clone() const
  {
    LinearElastic2DLaw::Pointer p_clone(new LinearElastic2DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LinearElastic2DLaw::~LinearElastic2DLaw()
  {
  }

  
  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool LinearElastic2DLaw::Has( const Variable<double>& rThisVariable )
  {
    return false;
  }

  bool LinearElastic2DLaw::Has( const Variable<Vector>& rThisVariable )
  {
    return false;
  }

  bool LinearElastic2DLaw::Has( const Variable<Matrix>& rThisVariable )
  {
    return false;
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& LinearElastic2DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    rValue = 0;

    if (rThisVariable==DETERMINANT_F)
      {
	rValue=mDetF0;
      }
  
    return rValue; 
  }

  Vector& LinearElastic2DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
  {

    if (rThisVariable==PK2_STRESS_TENSOR)
      {
	rValue =  mStressVector; 

	if(mStressMeasure!=StressMeasure_PK2)
	  KRATOS_WATCH(" this is not the PK2_STRESS_TENSOR ");
      }

    if (rThisVariable==CAUCHY_STRESS_TENSOR)
      {
	rValue = mStressVector;

	if(mStressMeasure!=StressMeasure_Cauchy)
	  KRATOS_WATCH(" this is not the CAUCHY_STRESS_TENSOR ");

      }

    /////**********///// transfer purpouses

    if (rThisVariable==CAUCHY_STRESS_VECTOR)
      {
	rValue = mStressVector;

	if(mStressMeasure!=StressMeasure_Cauchy)
	  KRATOS_WATCH(" this is not the CAUCHY_STRESS_VECTOR ");

      }

    if (rThisVariable==PK2_STRESS_VECTOR)
      {
	rValue = mStressVector;

	if(mStressMeasure!=StressMeasure_PK2)
	  KRATOS_WATCH(" this is not the PK2_STRESS_VECTOR ");

      }
				

    return( rValue );
  }

  Matrix& LinearElastic2DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
    if (rThisVariable==PK2_STRESS_TENSOR)
      {
	rValue = MathUtils<double>::StressVectorToTensor( mStressVector );

	if(mStressMeasure!=StressMeasure_PK2){
	  KRATOS_WATCH(" this is not the PK2_STRESS_TENSOR ");
	  KRATOS_WATCH(mStressVector);
	}

      }

    if (rThisVariable==CAUCHY_STRESS_TENSOR)
      {	        
	rValue = MathUtils<double>::StressVectorToTensor( mStressVector );

	if(mStressMeasure!=StressMeasure_Cauchy)
	  KRATOS_WATCH(" this is not the CAUCHY_STRESS_TENSOR ");

      }

    if (rThisVariable==DEFORMATION_GRADIENT) //returns the deformation gradient to compute the green lagrange strain
      {
	rValue = mDeformationGradientF0;
      }
    
    return( rValue );
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************


  void LinearElastic2DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

    if (rThisVariable==DETERMINANT_F)
      {
	mDetF0=rValue;
      }


  }

  void LinearElastic2DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
		
    if (rThisVariable==PK2_STRESS_TENSOR)
      {
	mStressVector=rValue;
      }

    if (rThisVariable==CAUCHY_STRESS_TENSOR)
      {
	mStressVector=rValue;
      }

    if (rThisVariable==CAUCHY_STRESS_VECTOR)
      {
	mStressVector=rValue;
      }
		
  }

  void LinearElastic2DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

    if (rThisVariable==DEFORMATION_GRADIENT)
      {
	mDeformationGradientF0=rValue;
      }


  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::InitializeMaterial( const Properties& props,
					      const GeometryType& geom,
					      const Vector& ShapeFunctionsValues )
  {
    unsigned int VoigtSize = 3;
    if(props[AXISYMMETRIC_LAW] == true)
      VoigtSize=4;
	  
    mStressVector = ZeroVector(VoigtSize);
    mDetF0 = 1;
    mDeformationGradientF0=identity_matrix<double>( VoigtSize-1 );
  }

  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::InitializeSolutionStep( const Properties& props,
						  const GeometryType& geom, //this is just to give the array of nodes
						  const Vector& ShapeFunctionsValues,
						  const ProcessInfo& CurrentProcessInfo)
  {
    mStressMeasure=StressMeasure_PK2;
  }
		
  //************************************************************************************
  //************************************************************************************

	
  void LinearElastic2DLaw::FinalizeSolutionStep( const Properties& props,
						const GeometryType& geom, //this is just to give the array of nodes
						const Vector& ShapeFunctionsValues,
						const ProcessInfo& CurrentProcessInfo)
  {
    mStressMeasure=StressMeasure_Cauchy;
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  LinearElastic2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
  {

    //-----------------------------//
   
    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    const Matrix& DeformationGradientF  = rValues.GetDeformationGradientF();
    Vector& StressVector                = rValues.GetStressVector();
    //const Vector& StrainVector          = rValues.GetStrainVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    const Properties&   props           = rValues.GetMaterialProperties();
    double detF                         = rValues.GetDeterminantF(detF); //Determinant of F from the reference configuration 
      
    //-----------------------------//

    //0.- Voigt size
    Flags &Options=rValues.GetOptions();
    mVoigtSize = 3;

    if(Options.Is(ConstitutiveLaw::AXISYMMETRIC))
      mVoigtSize = 4;

    //check if the given deformation gradient is the total or the incremental one for PK2
    if(Options.Is(ConstitutiveLaw::TOTAL_DEFORMATION_GRADIENT)){
      mDeformationGradientF0 = identity_matrix<double>( mVoigtSize-1 );
      mDetF0 = 1;
    }

    //1.- Lame constants
    const double& YoungModulus          = props[YOUNG_MODULUS];
    const double& PoissonCoefficient    = props[POISSON_RATIO];


    //2.-Right Cauchy-Green tensor C
    Matrix DeformationGradientF0 = prod(DeformationGradientF,mDeformationGradientF0);
     
    //3.-Right Cauchy Green
    Matrix RightCauchyGreen = prod(trans(DeformationGradientF0),DeformationGradientF0);

    //4.-Compute Green-Lagrange Strain  E= 0.5*(FT*F-1)
    Vector TotalStrainVector( mVoigtSize );
    TotalStrainVector[0] = 0.5 * ( RightCauchyGreen( 0, 0 ) - 1.00 );
    TotalStrainVector[1] = 0.5 * ( RightCauchyGreen( 1, 1 ) - 1.00 );
    TotalStrainVector[2] = RightCauchyGreen( 0, 1 );


    //7.-Incremental form
   
    if( Options.Is( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){
	  
      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

      CalculateStress( TotalStrainVector, ConstitutiveMatrix, StressVector );		

    }
    else if(  Options.IsNot( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){

      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

    }

   
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Strain "<<TotalStrainVector<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;
		
  }


  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
  {
    CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector=rValues.GetStressVector();
    const Matrix& F     =rValues.GetDeformationGradientF();
    const double& detF  =MathUtils<double>::Det(F);
    
    TransformStresses(StressVector,F,detF,StressMeasure_PK2,StressMeasure_PK1);
  }

  //************************************************************************************
  //************************************************************************************

  
  void LinearElastic2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
  {

    //-----------------------------//
   
    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    const Matrix& DeformationGradientF  = rValues.GetDeformationGradientF();
    Vector& StressVector                = rValues.GetStressVector();
    //const Vector& StrainVector          = rValues.GetStrainVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    const Properties&   props           = rValues.GetMaterialProperties();
    double detF                         = rValues.GetDeterminantF(detF); //Determinant of F from the reference configuration 
      
    //-----------------------------//

    //0.- Voigt size
    Flags &Options=rValues.GetOptions();
    mVoigtSize = 3;

    if(Options.Is(ConstitutiveLaw::AXISYMMETRIC))
      mVoigtSize = 4;

    //Set configuration where the the law is integrated:
    if(Options.Is(ConstitutiveLaw::INITIAL_CONFIGURATION)){
      mDeformationGradientF0 = identity_matrix<double>( mVoigtSize-1 );
      mDetF0 = 1;
    }

    //1.- Lame constants
    const double& YoungModulus          = props[YOUNG_MODULUS];
    const double& PoissonCoefficient    = props[POISSON_RATIO];


    //2.-Right Cauchy-Green tensor C
    Matrix DeformationGradientF0 = prod(DeformationGradientF,mDeformationGradientF0);

    //3.-Left Cauchy-Green tensor b
    Matrix LeftCauchyGreen = prod(DeformationGradientF0,trans(DeformationGradientF0));

    //4.-Invert the Left Cauchy-Green tensor b
    Matrix InverseLeftCauchyGreen ( mVoigtSize-1 , mVoigtSize-1 );
    double Trace_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, Trace_b);

    //5.-Compute almansi-strain  e= 0.5*(1-invbT*invb)
    Vector TotalStrainVector( mVoigtSize );
    TotalStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    TotalStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    TotalStrainVector[2] = InverseLeftCauchyGreen( 0, 1 );


    //7.-Incremental form
   
    if( Options.Is( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){
	  
      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

      CalculateStress( TotalStrainVector, ConstitutiveMatrix, StressVector );		

    }
    else if(  Options.IsNot( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){

      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

    }
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElastic2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
  {

    CalculateMaterialResponseKirchhoff (rValues);

    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    double detF                         = rValues.GetDeterminantF();
    double detF0                        = mDetF0 * detF;

    //Set to cauchy Stress:
    StressVector       /= detF0;
    ConstitutiveMatrix /= detF0;
  
    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;

  }


  //***********************************UPDATE*******************************************
  //************************************************************************************

  void LinearElastic2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
  {
	    
    CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =MathUtils<double>::Det(DeformationGradientF);

    //0.-Transform historical stresses to the current magnitude
    StressVector -= mStressVector;
    mStressVector+= StressVector;

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_Cauchy); //total Cauchy Stress

    //2.-Returns the stress increment to the updated configuration
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //3.-Update Internal Variables
    mDeformationGradientF0  = prod(DeformationGradientF,mDeformationGradientF0);
    mDetF0 *= detF;
  }

  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
  {
	    
    CalculateMaterialResponsePK1 (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =MathUtils<double>::Det(DeformationGradientF);


    //0.-Transform historical stresses to the current magnitude
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_PK1);

    StressVector -= mStressVector;
    mStressVector+= StressVector;

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_PK1,StressMeasure_Cauchy); //total Cauchy Stress

    //2.-Returns the stress increment to the updated configuration
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_PK1,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //3.-Update Internal Variables
    mDeformationGradientF0  = prod(DeformationGradientF,mDeformationGradientF0);
    mDetF0 *= detF;
  }

  //************************************************************************************
  //************************************************************************************

  
  void LinearElastic2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
  {		
   
    CalculateMaterialResponseKirchhoff (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =MathUtils<double>::Det(DeformationGradientF);

    //0.-Transform historical stresses to the current magnitude
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_Kirchhoff);

    StressVector -= mStressVector;
    mStressVector+= StressVector;

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_Kirchhoff,StressMeasure_Cauchy); //total Cauchy Stress

    //2.-Returns the stress increment to the updated configuration
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_Kirchhoff,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //3.-Update Internal Variables
    mDeformationGradientF0  = prod(DeformationGradientF,mDeformationGradientF0);
    mDetF0 *= detF;
  }


  //************************************************************************************
  //************************************************************************************

  void LinearElastic2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
  {

   
    CalculateMaterialResponseCauchy (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =MathUtils<double>::Det(DeformationGradientF);

    //0.-Transform historical stresses to the current magnitude
    TransformStresses(mStressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_Cauchy);

    StressVector -= mStressVector;
    mStressVector+= StressVector;

    //3.-Update Internal Variables
    mDeformationGradientF0  = prod(DeformationGradientF,mDeformationGradientF0);
    mDetF0 *= detF;

  }

  

  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void LinearElastic2DLaw::CalculateStress( const Vector & rStrainVector,
					   const Matrix & rConstitutiveMatrix,
					   Vector& rStressVector )
  {
    
    //1.-Calculate Stress Increment
    CalculateStressIncrement(rStrainVector,rConstitutiveMatrix,rStressVector);

    //2.-Add the hitorical 2nd Piola Kirchooff stress:
    for(unsigned int i = 0; i<mStressVector.size(); i++){
      rStressVector[i]+=mStressVector[i];
    }
    
  }

 

  //***********************COMPUTE STRESS INCREMENT ************************************
  //************************************************************************************

  void LinearElastic2DLaw::CalculateStressIncrement(const Vector & rStrainVector,
						    const Matrix & rConstitutiveMatrix,
						    Vector& rStressVector)
  {
	  
    //1.-2nd Piola Kirchhoff StressVector increment
    rStressVector = prod(rConstitutiveMatrix,rStrainVector);

    //2.-Increment of Stresses:
    for(unsigned int i = 0; i<rStressVector.size(); i++){
      rStressVector[i]-=mStressVector[i];
    }
    
  }



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************


  void LinearElastic2DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, 
							const double &rYoungModulus, 
							const double &rPoissonCoefficient )
  {
    rConstitutiveMatrix.clear();

    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1-2*rPoissonCoefficient)/(2*(1.0-rPoissonCoefficient));

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );
  }


  

 
  //*******************************METHOD FROM BASE CLASS******************************
  //************************************************************************************
    
  void LinearElastic2DLaw::CalculateCauchyStresses(Vector& rCauchy_StressVector,
						  const Matrix& rF,
						  const Vector& rPK2_StressVector,
						  const Vector& rGreenLagrangeStrainVector )
  {

  }



  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  bool LinearElastic2DLaw::CheckParameters(Parameters& rValues)
  {
    return rValues.CheckAllParameters();
  }



  int LinearElastic2DLaw::Check(const Properties& rProperties, 
			       const GeometryType& rGeometry, 
			       const ProcessInfo& rCurrentProcessInfo)
  {

    if(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS]<= 0.00)
      KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = rProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true) 
      KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

	  
    if(DENSITY.Key() == 0 || rProperties[DENSITY]<0.00)
      KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");

	  	    
    return 0;
	    
  }
    
} // Namespace Kratos
