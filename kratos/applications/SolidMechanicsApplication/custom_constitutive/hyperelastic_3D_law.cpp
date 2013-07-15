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

  HyperElastic3DLaw::HyperElastic3DLaw()
  : ConstitutiveLaw()
  {

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElastic3DLaw::HyperElastic3DLaw(const HyperElastic3DLaw& rOther)
  : ConstitutiveLaw()
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElastic3DLaw::Clone() const
  {
    HyperElastic3DLaw::Pointer p_clone(new HyperElastic3DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElastic3DLaw::~HyperElastic3DLaw()
  {
  }

  
  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool HyperElastic3DLaw::Has( const Variable<double>& rThisVariable )
  {
    return false;
  }

  bool HyperElastic3DLaw::Has( const Variable<Vector>& rThisVariable )
  {
    return false;
  }

  bool HyperElastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
  {
    return false;
  }


  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& HyperElastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {

    return( rValue ); 
  }

  Vector& HyperElastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
  {

    return( rValue );
  }

  Matrix& HyperElastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
 
    return( rValue );
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************


  void HyperElastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

  }

  void HyperElastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {
		
		
  }

  void HyperElastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void HyperElastic3DLaw::InitializeMaterial( const Properties& props,
					      const GeometryType& geom,
					      const Vector& ShapeFunctionsValues )
  {
   
  }

  //************************************************************************************
  //************************************************************************************


  void HyperElastic3DLaw::InitializeSolutionStep( const Properties& props,
						  const GeometryType& geom, //this is just to give the array of nodes
						  const Vector& ShapeFunctionsValues,
						  const ProcessInfo& CurrentProcessInfo)
  {
    
  }
		
  //************************************************************************************
  //************************************************************************************

	
  void HyperElastic3DLaw::FinalizeSolutionStep( const Properties& props,
						const GeometryType& geom, //this is just to give the array of nodes
						const Vector& ShapeFunctionsValues,
						const ProcessInfo& CurrentProcessInfo)
  {
    
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  HyperElastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
  {

    //-----------------------------//
   
    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& detF                    = rValues.GetDeterminantF(); 

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& detF0                         = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Dimension
    const unsigned int dimension =  DeformationGradientF.size1();

    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);

    //3.-Determinant of the Total Deformation Gradient
    detF0 *= detF;
   
    //4.-Right Cauchy Green
    Matrix RightCauchyGreen = prod(trans(F0),F0);

    //5.-Inverse of the Right Cauchy-Green tensor C:
    double Trace_C=0;
    Matrix InverseRightCauchyGreen ( dimension , dimension );
    MathUtils<double>::InvertMatrix( RightCauchyGreen, InverseRightCauchyGreen, Trace_C);

    //6.-Green-Lagrange Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);
      }

    //7.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( dimension );


    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
		  
 	this->CalculateStress( InverseRightCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_PK2, StressVector );
    }
   
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

 	this->CalculateConstitutiveMatrix ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, ConstitutiveMatrix );
    }
    
    //OPTION 2:
    if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) ){

      if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
	//Left Cauchy-Green tensor b
	Matrix LeftCauchyGreen = prod(F0,trans(F0));

	this->CalculateStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_Kirchhoff, StressVector );

	TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_Kirchhoff,StressMeasure_PK2);  //2nd PK Stress in the last known configuration
      }

      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
     
	Matrix invF ( dimension , dimension );
	double DetInvF=0;
	MathUtils<double>::InvertMatrix( DeformationGradientF, invF, DetInvF);
		
	this->CalculateConstitutiveMatrix (  IdentityMatrix, invF, detF0, LameLambda, LameMu, ConstitutiveMatrix );

	ConstitutiveMatrix *= detF;
      }

    }
  

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;
		
  }


  //************************************************************************************
  //************************************************************************************


  void HyperElastic3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
  {
    this->CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector=rValues.GetStressVector();
    const Matrix& F     =rValues.GetDeformationGradientF();
    const double& detF  =rValues.GetDeterminantF();
    
    TransformStresses(StressVector,F,detF,StressMeasure_PK2,StressMeasure_PK1);
  }

  //************************************************************************************
  //************************************************************************************

  
  void HyperElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
  {

    //-----------------------------//
   
    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   detF                  = rValues.GetDeterminantF(); 

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& detF0                         = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Dimension
    const unsigned int dimension =  DeformationGradientF.size1();

    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);

    //3.-Determinant of the Total Deformation Gradient
    detF0 *= detF;
        
    //4.-Left Cauchy-Green tensor b
    Matrix LeftCauchyGreen = prod(F0,trans(F0));

    //6.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	// e= 0.5*(1-invbT*invb)   
	this->CalculateAlmansiStrain(LeftCauchyGreen,StrainVector);
      }
 
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( dimension );

 		
    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
      this->CalculateStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_Kirchhoff, StressVector );
   
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateConstitutiveMatrix ( IdentityMatrix, detF0, LameLambda, LameMu, ConstitutiveMatrix );
         
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
  {

    this->CalculateMaterialResponseKirchhoff (rValues);

    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    double& detF0                       = rValues.GetDeterminantF0();;

    //Set to cauchy Stress:
    StressVector       /= detF0;
    ConstitutiveMatrix /= detF0;
  
    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;

  }


  //***********************************UPDATE*******************************************
  //************************************************************************************

  void HyperElastic3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
  {
	    
    this->CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& detF0                          =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_PK2,StressMeasure_Cauchy);  //Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    detF0                  = MathUtils<double>::Det(DeformationGradientF0);
  }

  //************************************************************************************
  //************************************************************************************


  void HyperElastic3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
  {
	    
    this->CalculateMaterialResponsePK1 (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& detF0                          =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step  
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_PK1,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    detF0                  = MathUtils<double>::Det(DeformationGradientF0);
  }

  //************************************************************************************
  //************************************************************************************

  
  void HyperElastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
  {		
   
    this->CalculateMaterialResponseKirchhoff (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& detF0                          =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& detF                     =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step  
    TransformStresses(StressVector,DeformationGradientF,detF,StressMeasure_Kirchhoff,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    detF0                  = MathUtils<double>::Det(DeformationGradientF0);
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
  {

   
    this->CalculateMaterialResponseCauchy (rValues);
  
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& detF0                          =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    detF0                  = MathUtils<double>::Det(DeformationGradientF0);


  }



  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
							Vector& rStrainVector )
  {

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );   
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

  }



  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						  Vector& rStrainVector )
  {

    // e= 0.5*(1-invbT*invb)  
 
   //Calculating the inverse of the jacobian 
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz 
  }

  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateStress( const Matrix & rMatrixIC,
					   const Matrix & rIdentityMatrix,
					   const double & rdetF0,
					   const double & rLameLambda, 
					   const double & rLameMu, 
					   StressMeasure rStressMeasure,
					   Vector& rStressVector )
  {

    //1.- Temporary and selected law
    Matrix StressMatrix   ( rMatrixIC.size1() , rMatrixIC.size1() );

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


  }


  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    //diagonal
    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,1,1);
    //C3333
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,2,2);
    //C1212
    rConstitutiveMatrix( 3, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,1);
    //C2323
    rConstitutiveMatrix( 4, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,1,2);
    //C1313
    rConstitutiveMatrix( 5, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,2,0,2);

    //row 1
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,1,1); 
    //C1133
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,2,2);
    //C1112
    rConstitutiveMatrix( 0, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,1);
    //C1123
    rConstitutiveMatrix( 0, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,1,2);
    //C1113
    rConstitutiveMatrix( 0, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,0,0,2);
    
    //row 2
    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,0);
    //C2233
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,2,2);
    //C2212
    rConstitutiveMatrix( 1, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,1);
    //C2223
    rConstitutiveMatrix( 1, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,1,2);
    //C2213
    rConstitutiveMatrix( 1, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,1,0,2);

    //row 3
    //C3311
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,0,0);
    //C3322
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,1,1);
    //C3312
    rConstitutiveMatrix( 2, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,0,1);
    //C3323
    rConstitutiveMatrix( 2, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,1,2);
    //C3313
    rConstitutiveMatrix( 2, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,2,2,0,2);

    //row 4
    //C1211
    rConstitutiveMatrix( 3, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    //C1222
    rConstitutiveMatrix( 3, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1233
    rConstitutiveMatrix( 3, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,2,2);
    //C1223
    rConstitutiveMatrix( 3, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,2);
    //C1213
    rConstitutiveMatrix( 3, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,2);

    //row 5
    //C2311
    rConstitutiveMatrix( 4, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,0,0);
    //C2322
    rConstitutiveMatrix( 4, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,1,1);
    //C2333
    rConstitutiveMatrix( 4, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,2,2);
    //C2312
    rConstitutiveMatrix( 4, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,0,1);
    //C2313
    rConstitutiveMatrix( 4, 5 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,1,2,0,2);

    //row 6
    //C1311
    rConstitutiveMatrix( 5, 0 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    //C1322
    rConstitutiveMatrix( 5, 1 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1333
    rConstitutiveMatrix( 5, 2 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,2,2);
    //C1312
    rConstitutiveMatrix( 5, 3 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,0,1);
    //C1323
    rConstitutiveMatrix( 5, 4 )=ConstitutiveComponent(rMatrixIC,rdetF0,rLameLambda,rLameMu,0,1,1,2);

	  	
  }


  //**************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX PULL-BACK*********************
  //************************************************************************************

  void HyperElastic3DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const Matrix &rinvF,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rConstitutiveMatrix)
  {

    rConstitutiveMatrix.clear();
		
    //diagonal
    //C1111
    rConstitutiveMatrix( 0, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,0);
    //C2222
    rConstitutiveMatrix( 1, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,1,1);
    //C3333
    rConstitutiveMatrix( 2, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,2,2);
    //C1212
    rConstitutiveMatrix( 3, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,1);
    //C2323
    rConstitutiveMatrix( 4, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,1,2);
    //C1313
    rConstitutiveMatrix( 5, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,2,0,2);

    //row 1
    //C1122
    rConstitutiveMatrix( 0, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,1,1); 
    //C1133
    rConstitutiveMatrix( 0, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,2,2);
    //C1112
    rConstitutiveMatrix( 0, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,1);
    //C1123
    rConstitutiveMatrix( 0, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,1,2);
    //C1113
    rConstitutiveMatrix( 0, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,0,0,2);
    
    //row 2
    //C2211
    rConstitutiveMatrix( 1, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,0);
    //C2233
    rConstitutiveMatrix( 1, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,2,2);
    //C2212
    rConstitutiveMatrix( 1, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,1);
    //C2223
    rConstitutiveMatrix( 1, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,1,2);
    //C2213
    rConstitutiveMatrix( 1, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,1,0,2);

    //row 3
    //C3311
    rConstitutiveMatrix( 2, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,0,0);
    //C3322
    rConstitutiveMatrix( 2, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,1,1);
    //C3312
    rConstitutiveMatrix( 2, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,0,1);
    //C3323
    rConstitutiveMatrix( 2, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,1,2);
    //C3313
    rConstitutiveMatrix( 2, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,2,2,0,2);

    //row 4
    //C1211
    rConstitutiveMatrix( 3, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    //C1222
    rConstitutiveMatrix( 3, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1233
    rConstitutiveMatrix( 3, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,2,2);
    //C1223
    rConstitutiveMatrix( 3, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,2);
    //C1213
    rConstitutiveMatrix( 3, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,2);

    //row 5
    //C2311
    rConstitutiveMatrix( 4, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,0,0);
    //C2322
    rConstitutiveMatrix( 4, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,1,1);
    //C2333
    rConstitutiveMatrix( 4, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,2,2);
    //C2312
    rConstitutiveMatrix( 4, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,0,1);
    //C2313
    rConstitutiveMatrix( 4, 5 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,1,2,0,2);

    //row 6
    //C1311
    rConstitutiveMatrix( 5, 0 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,0);
    //C1322
    rConstitutiveMatrix( 5, 1 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,1);
    //C1333
    rConstitutiveMatrix( 5, 2 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,2,2);
    //C1312
    rConstitutiveMatrix( 5, 3 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,0,1);
    //C1323
    rConstitutiveMatrix( 5, 4 )=ConstitutiveComponent(rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,0,1,1,2);

    	  
  }




  
  //***********************CONSTITUTIVE MATRIX COMPONENTS*******************************
  //************************************************************************************


  double HyperElastic3DLaw::ConstitutiveComponent(const Matrix &rMatrixIC, 
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


  double HyperElastic3DLaw::ConstitutiveComponent(const Matrix &rMatrixIC, 
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
    int dim = rinvF.size1();
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


  //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
  //************************************************************************************

  bool HyperElastic3DLaw::CheckParameters(Parameters& rValues)
  {
    return rValues.CheckAllParameters();
  }



  int HyperElastic3DLaw::Check(const Properties& rProperties, 
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
