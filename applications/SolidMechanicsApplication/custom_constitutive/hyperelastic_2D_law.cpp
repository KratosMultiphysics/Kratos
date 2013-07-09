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
#include "custom_constitutive/hyperelastic_2D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElastic2DLaw::HyperElastic2DLaw()
  : ConstitutiveLaw()
  {
    //pointer to constitutive law in properties :: serializer has nan in not set internal variables
    mDetF0 = 0;

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElastic2DLaw::HyperElastic2DLaw(const HyperElastic2DLaw& rOther)
  : ConstitutiveLaw()
    ,mStressVector(rOther.mStressVector)
    ,mStressMeasure(rOther.mStressMeasure)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDetF0(rOther.mDetF0)
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElastic2DLaw::Clone() const
  {
    HyperElastic2DLaw::Pointer p_clone(new HyperElastic2DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElastic2DLaw::~HyperElastic2DLaw()
  {
  }

  
  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool HyperElastic2DLaw::Has( const Variable<double>& rThisVariable )
  {
    return false;
  }

  bool HyperElastic2DLaw::Has( const Variable<Vector>& rThisVariable )
  {
    return false;
  }

  bool HyperElastic2DLaw::Has( const Variable<Matrix>& rThisVariable )
  {
    return false;
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& HyperElastic2DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    rValue = 0;

    if (rThisVariable==DETERMINANT_F)
      {
	rValue=mDetF0;
      }
  
    return rValue; 
  }

  Vector& HyperElastic2DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
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

  Matrix& HyperElastic2DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
    if (rThisVariable==PK2_STRESS_TENSOR)
      {
	rValue = MathUtils<double>::StressVectorToTensor( mStressVector );

	if(mStressMeasure!=StressMeasure_PK2)
	  KRATOS_WATCH(" this is not the PK2_STRESS_TENSOR ");

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


  void HyperElastic2DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

    if (rThisVariable==DETERMINANT_F)
      {
	mDetF0=rValue;
      }


  }

  void HyperElastic2DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
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

  void HyperElastic2DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
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


  void HyperElastic2DLaw::InitializeMaterial( const Properties& props,
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


  void HyperElastic2DLaw::InitializeSolutionStep( const Properties& props,
						  const GeometryType& geom, //this is just to give the array of nodes
						  const Vector& ShapeFunctionsValues,
						  const ProcessInfo& CurrentProcessInfo)
  {
    mStressMeasure=StressMeasure_PK2;
  }
		
  //************************************************************************************
  //************************************************************************************

	
  void HyperElastic2DLaw::FinalizeSolutionStep( const Properties& props,
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


  void  HyperElastic2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
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

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Right Cauchy-Green tensor C
    Matrix DeformationGradientF0 = prod(DeformationGradientF,mDeformationGradientF0);

    //3.-Determinant Deformation Gradient
    double detF0 = mDetF0 * detF;
   
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( mVoigtSize-1 );

    
    //5.-Right Cauchy Green
    Matrix RightCauchyGreen = prod(trans(DeformationGradientF0),DeformationGradientF0);

    //6.-Inverse of the Right Cauchy-Green tensor C:
    Matrix InverseRightCauchyGreen ( mVoigtSize-1 , mVoigtSize-1 );
    double Trace_C=0;
    MathUtils<double>::InvertMatrix( RightCauchyGreen, InverseRightCauchyGreen, Trace_C);

    //7.-Incremental form

    //OPTION 1:
    if( Options.Is( COMPUTE_STRESS ) ){
		  
      if( Options.Is( LAST_KNOWN_CONFIGURATION ) ){
	//Left Cauchy-Green tensor b
	Matrix LeftCauchyGreen = prod(DeformationGradientF0,trans(DeformationGradientF0));
	CalculateStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_Kirchhoff, StressVector );

	TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_Kirchhoff,StressMeasure_PK2);  //2nd PK Stress
      }
      else{

	CalculateStress( InverseRightCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_PK2, StressVector );
      }
    }
   
    if( Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){

      if( Options.Is( LAST_KNOWN_CONFIGURATION ) ){
	Matrix invF ( mVoigtSize-1, mVoigtSize-1 );
	double DetInvF=0;
	MathUtils<double>::InvertMatrix( DeformationGradientF, invF, DetInvF);
		
	CalculateConstitutiveMatrix (  IdentityMatrix, invF, detF0, LameLambda, LameMu, ConstitutiveMatrix );

	ConstitutiveMatrix *= detF;
      }
      else{
	CalculateConstitutiveMatrix ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, ConstitutiveMatrix );
      }

    }
    
    //OPTION 2:
    // if( Options.Is( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) )
    // {
    //   CalculateConstitutiveMatrix ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, ConstitutiveMatrix );
    //   CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );		
    // }
    
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;
		
  }


  //************************************************************************************
  //************************************************************************************


  void HyperElastic2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
  {
    CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector=rValues.GetStressVector();
    const Matrix& F     =rValues.GetDeformationGradientF();
    const double& detF  =MathUtils<double>::Det(F);
    
    TransformStresses(StressVector,F,detF,StressMeasure_PK2,StressMeasure_PK1);
  }

  //************************************************************************************
  //************************************************************************************

  
  void HyperElastic2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
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
      mVoigtSize = 2;

    //Set configuration where the the law is integrated:
    if(Options.Is(ConstitutiveLaw::INITIAL_CONFIGURATION)){
      mDeformationGradientF0 = identity_matrix<double>( mVoigtSize-1 );
      mDetF0 = 1;
    }

    //1.- Lame constants
    const double& YoungModulus          = props[YOUNG_MODULUS];
    const double& PoissonCoefficient    = props[POISSON_RATIO];

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Right Cauchy-Green tensor C
    Matrix DeformationGradientF0 = prod(DeformationGradientF,mDeformationGradientF0);

    //3.-Determinant Deformation Gradient
    double detF0 = mDetF0 * detF;
   
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( mVoigtSize-1 );

    //5.-Left Cauchy-Green tensor b
    Matrix LeftCauchyGreen = prod(DeformationGradientF0,trans(DeformationGradientF0));
		
    //OPTION 1:
    if( Options.Is( COMPUTE_STRESS ) )
      CalculateStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_Kirchhoff, StressVector );
   
    if( Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) )
      CalculateConstitutiveMatrix ( IdentityMatrix, detF0, LameLambda, LameMu, ConstitutiveMatrix );
    
    //OPTION 2:
    // if( Options.Is( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) )
    // {
    //   CalculateConstitutiveMatrix ( IdentityMatrix, detF0, LameLambda, LameMu, ConstitutiveMatrix );
    //   CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );		
    // }
      
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
  {

    CalculateMaterialResponseKirchhoff (rValues);

    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    double detF                         = rValues.GetDeterminantF();
    double detF0                        = mDetF0 * detF;

    //Set to cauchy Stress:
    StressVector       /= detF0;
    ConstitutiveMatrix /= detF0;
  
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;

  }


  //***********************************UPDATE*******************************************
  //************************************************************************************

  void HyperElastic2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
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


  void HyperElastic2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
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

  
  void HyperElastic2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
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

  void HyperElastic2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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

  void HyperElastic2DLaw::CalculateStress( const Matrix & rMatrixIC,
					   const Matrix & rIdentityMatrix,
					   const double & rdetF0,
					   const double & rLameLambda, 
					   const double & rLameMu, 
					   StressMeasure rStressMeasure,
					   Vector& rStressVector )
  {

    //1.-Calculate Stress Increment
    CalculateStressIncrement(rMatrixIC,rIdentityMatrix,rdetF0,rLameLambda,rLameMu,rStressMeasure,rStressVector);

    //2.-Add the hitorical 2nd Piola Kirchooff stress:
    for(unsigned int i = 0; i<mStressVector.size(); i++){
      rStressVector[i]+=mStressVector[i];
    }
    
    //std::cout<<" StressVector "<<rStressVector<<std::endl;
    //std::cout<<" rMatrixIC "<<rMatrixIC<<std::endl;
  }



  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void HyperElastic2DLaw::CalculateStress( const Vector & rStrainVector,
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

  void HyperElastic2DLaw::CalculateStressIncrement( const Matrix & rMatrixIC,
						    const Matrix & rIdentityMatrix,
						    const double & rdetF0,
						    const double & rLameLambda, 
						    const double & rLameMu, 
						    StressMeasure rStressMeasure,
						    Vector& rStressVector )
  {
    //1.- Temporary and selected law
    Matrix StressMatrix   ( mVoigtSize-1 , mVoigtSize-1 );

    double auxiliar = (std::log(rdetF0)); //(ln(J))
    //double auxiliar = 0.5*(rdetF0*rdetF0-1); //(J²-1)/2

    if(rStressMeasure == StressMeasure_PK2){    
    
      //rMatrixIC is InverseRightCauchyGreen

      //2.-2nd Piola Kirchhoff Stress Matrix
      StressMatrix  = rLameLambda*auxiliar*rMatrixIC;
      StressMatrix += rLameMu*(rIdentityMatrix-rMatrixIC);
    }

    if(rStressMeasure == StressMeasure_Kirchhoff){

      //rMatrixIC is LeftCauchyGreen
		  
      //2.-Kirchhoff Stress Matrix 
      StressMatrix  = rLameLambda*auxiliar*rIdentityMatrix;
      StressMatrix += rLameMu*(rMatrixIC-rIdentityMatrix);
    }

    //3.-Add the hitorical 2nd Piola Kirchooff stress:
    rStressVector[0]=StressMatrix( 0 , 0 )-mStressVector[0];
    rStressVector[1]=StressMatrix( 1 , 1 )-mStressVector[1];

    if(rStressVector.size() == 4){
      rStressVector[2]=StressMatrix( 2 , 2 )-mStressVector[2];
      rStressVector[3]=StressMatrix( 0 , 1 )-mStressVector[3];
    }
    else{
      rStressVector[2]=StressMatrix( 0 , 1 )-mStressVector[2];
    }

		
  }

  

  //***********************COMPUTE STRESS INCREMENT ************************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateStressIncrement(const Vector & rStrainVector,
						   const Matrix & rConstitutiveMatrix,
						   Vector& rStressVector)
  {
	  
    //1.-2nd Piola Kirchhoff StressVector increment
    rStressVector = prod(rConstitutiveMatrix,rStrainVector);
  
  }



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rResult)
  {

    if( mVoigtSize == 4)
      CalculateAxisymNeoHookeanMatrix(rResult,rMatrixIC,rdetF0,rLameLambda,rLameMu);
    else
      CalculateNeoHookeanMatrix(rResult,rMatrixIC,rdetF0,rLameLambda,rLameMu);
	  
    // double YoungModulus          = rLameMu*(3*rLameLambda+2*rLameMu)/(rLameLambda+rLameMu);
    // double PoissonCoefficient    = 0.5*rLameLambda/(rLameLambda+rLameMu);
	  
    // CalculateLinearElasticMatrix(rResult,YoungModulus,PoissonCoefficient);

    // std::cout<<" ConstitutiveTensor "<<rResult<<std::endl;
	
  }


  //**************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX PULL-BACK*********************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const Matrix &rinvF,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rResult)
  {
	  
    if( mVoigtSize == 4)
      CalculateAxisymNeoHookeanMatrix(rResult,rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu);
    else
      CalculateNeoHookeanMatrix(rResult,rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu);
	  
    //std::cout<<" ConstitutiveTensor "<<rResult<<std::endl;

  }



  //*******************************NEOHOOKEAN MATRIX************************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateNeoHookeanMatrix( Matrix& rConstitutiveMatrix, 
						     const Matrix &rMatrixIC, 
						     const double &rdetF0, 
						     const double &rLameLambda, 
						     const double &rLameMu )
  {
        
    rConstitutiveMatrix.clear();
		
    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,1,1);	       		//C1212
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,1);
   
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,1,1);
    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,0);

    //C1112
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,1);
    
    //C2212
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,1);

    //C1211
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    
    //C1222
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,1);
        
  }



  //*******************************NEOHOOKEAN MATRIX PULL-BACK**************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateNeoHookeanMatrix( Matrix& rConstitutiveMatrix, 
						     const Matrix &rMatrixIC, 
						     const Matrix &rinvF,
						     const double &rdetF0, 
						     const double &rLameLambda, 
						     const double &rLameMu )
  {
        
    rConstitutiveMatrix.clear();

    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,1,1);
    //C1212
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,1);
   
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,1,1);
    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,0);

    //C1112
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,1);
    
    //C2212
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,1);

    //C1211
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    
    //C1222
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,1);
        
  }


  //*******************************NEOHOOKEAN MATRIX AXISYM*****************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateAxisymNeoHookeanMatrix( Matrix& rConstitutiveMatrix, 
							   const Matrix &rMatrixIC, 
							   const double &rdetF0, 
							   const double &rLameLambda, 
							   const double &rLameMu )
  {
        
    rConstitutiveMatrix.clear();
		
    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,1,1);	
    //C3333
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,2,2);	

    //C1212
    rConstitutiveMatrix( 3, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,1);
   
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,1,1);
    //C1133
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,2,2);

    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,0);
    //C2233
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,2,2);

    //C3311
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,0,0);
    //C3322
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,1,1);

    //C1112
    rConstitutiveMatrix( 0, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,1); 
    //C2212
    rConstitutiveMatrix( 1, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,1);
    //C3312
    rConstitutiveMatrix( 2, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,0,1);

    //C1211
    rConstitutiveMatrix( 3, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,0);   
    //C1222
    rConstitutiveMatrix( 3, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1233
    rConstitutiveMatrix( 3, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,2,2);

        
  }



  //***************************NEOHOOKEAN MATRIX AXISYM PULL-BACK***********************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateAxisymNeoHookeanMatrix( Matrix& rConstitutiveMatrix, 
							   const Matrix &rMatrixIC, 
							   const Matrix &rinvF,
							   const double &rdetF0, 
							   const double &rLameLambda, 
							   const double &rLameMu )
  {
        
    rConstitutiveMatrix.clear();

    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,1,1);	
    //C3333
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,2,2);	

    //C1212
    rConstitutiveMatrix( 3, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,1);
   
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,1,1);
    //C1133
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,2,2);

    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,0);
    //C2233
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,2,2);

    //C3311
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,0,0);
    //C3322
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,1,1);

    //C1112
    rConstitutiveMatrix( 0, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,1); 
    //C2212
    rConstitutiveMatrix( 1, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,1);
    //C3312
    rConstitutiveMatrix( 2, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,0,1);

    //C1211
    rConstitutiveMatrix( 3, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,0);   		//C1222
    rConstitutiveMatrix( 3, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1233
    rConstitutiveMatrix( 3, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,2,2);


  }



  void HyperElastic2DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, 
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


  
  //***********************CONSTITUTIVE MATRIX COMPONENTS*******************************
  //************************************************************************************


  double HyperElastic2DLaw::ConstitutiveComponent(const Matrix &rMatrixIC, 
						  const double &rdetF0, 
						  const double &rLameLambda, 
						  const double &rLameMu, 
						  int a, int b, int c, int d)
  {
	  
    //(J²-1)/2
    //double auxiliar1 =  rdetF0*rdetF0;
    //double auxiliar2 =  (rdetF0*rdetF0-1);

    //(ln(J))
    double auxiliar1 =  1.0/rdetF0;
    double auxiliar2 =  (2.0*std::log(rdetF0));


    //1.Elastic constitutive tensor component:
    double Cabcd=(rLameLambda*auxiliar1*rMatrixIC(a,b)*rMatrixIC(c,d));
    Cabcd+=((2*rLameMu-rLameLambda*auxiliar2)*0.5*(rMatrixIC(a,c)*rMatrixIC(b,d)+rMatrixIC(a,d)*rMatrixIC(b,c)));

    return Cabcd;
  }



  //***********************CONSTITUTIVE MATRIX COMPONENTS PULL-BACK*********************
  //************************************************************************************


  double HyperElastic2DLaw::ConstitutiveComponent(const Matrix &rMatrixIC, 
						  const Matrix &rinvF,
						  const double &rdetF0, 
						  const double &rLameLambda, 
						  const double &rLameMu, 
						  int a, int b, int c, int d)

  {

    double component=0;

    Matrix invF = identity_matrix<double> (3);

    if(rinvF.size1() == 2)
      {	
	for(int r=0; r<2; r++)
	  {   
	    for(int s=0; s<2; s++)
	      {
		invF(r,s) = rinvF(r,s);
	      }
	  }
      }
    else{

      invF=rinvF;
    }
    
    //std::cout<<" invF "<<invF<<std::endl;
    int dim=mVoigtSize-1;
    //Cabcd
    for(int j=0; j<dim; j++)
      {	  
	for(int l=0; l<dim; l++)
	  {	   	    
	    for(int k=0; k<dim; k++)
	      {
		for(int i=0; i<dim; i++)
		  {
		    //Cijkl
		    component +=invF(a,i)*invF(b,j)*invF(c,k)*invF(d,l)*ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,i,j,k,l);
		  }
	      }
	  }
      }

    return component;

  }


 
  //*******************************METHOD FROM BASE CLASS******************************
  //************************************************************************************
    
  void HyperElastic2DLaw::CalculateCauchyStresses(Vector& rCauchy_StressVector,
						  const Matrix& rF,
						  const Vector& rPK2_StressVector,
						  const Vector& rGreenLagrangeStrainVector )
  {

  }



  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  bool HyperElastic2DLaw::CheckParameters(Parameters& rValues)
  {
    return rValues.CheckAllParameters();
  }



  int HyperElastic2DLaw::Check(const Properties& rProperties, 
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
