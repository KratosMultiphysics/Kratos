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
#include "custom_constitutive/hyperelastic_3D_law.hpp"

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

  Matrix& HyperElastic3DLaw::DeformationGradient3D (Matrix & Matrix2D)
  {
    //Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal

    if (Matrix2D.size1() == 2 && Matrix2D.size2() == 2){

      Matrix2D.resize( 3, 3, true);

      Matrix2D( 0 , 2 ) = 0.0;
      Matrix2D( 1 , 2 ) = 0.0;

      Matrix2D( 2 , 0 ) = 0.0;
      Matrix2D( 2 , 1 ) = 0.0;

      Matrix2D( 2 , 2 ) = 1.0;

    }
    else if(Matrix2D.size1() != 3 && Matrix2D.size2() != 3){

      KRATOS_ERROR(std::invalid_argument,"Passed Matrix dimensions in DeformtationGradient3D not correct ","");

    }

    return Matrix2D;
    
  }



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
    const double& DeterminantF            = rValues.GetDeterminantF(); 

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Dimension
    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);

    F0 = DeformationGradient3D( F0 );

    //3.-Determinant of the Total Deformation Gradient
    double detF0 = DeterminantF0 * DeterminantF;
   
    //4.-Right Cauchy Green
    Matrix RightCauchyGreen = prod(trans(F0),F0);

    //5.-Inverse of the Right Cauchy-Green tensor C:
    double Trace_C=0;
    Matrix InverseRightCauchyGreen ( 3 , 3 );
    MathUtils<double>::InvertMatrix( RightCauchyGreen, InverseRightCauchyGreen, Trace_C);

    //6.-Green-Lagrange Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);
      }

    //7.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( 3 );


    //OPTION 1: ( initial configuration )
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
		  
 	this->CalculateStress( InverseRightCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_PK2, StressVector );
    }
   
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

 	this->CalculateConstitutiveMatrix ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, ConstitutiveMatrix );
    }
    
    //OPTION 2: ( last known configuration : updated lagrangian approach only )
    if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) ){

      if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){
	//Left Cauchy-Green tensor b
	Matrix LeftCauchyGreen = prod(F0,trans(F0));

	this->CalculateStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, LameMu, StressMeasure_Kirchhoff, StressVector );

	TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_Kirchhoff,StressMeasure_PK2); //2nd PK Stress in the last known configuration
      }

      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
     
	Matrix invF ( 3, 3 );
	double DetInvF=0;
	MathUtils<double>::InvertMatrix( DeformationGradientF, invF, DetInvF);
		
	this->CalculateConstitutiveMatrix (  IdentityMatrix, invF, detF0, LameLambda, LameMu, ConstitutiveMatrix );

	ConstitutiveMatrix *= DeterminantF;
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
  
    Vector& StressVector = rValues.GetStressVector();
    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();
    
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_PK1);
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
    const double&   DeterminantF          = rValues.GetDeterminantF(); 

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Dimension
    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    double LameLambda = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu     =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);
    F0 = DeformationGradient3D( F0 );

    //3.-Determinant of the Total Deformation Gradient
    double detF0 = DeterminantF0 * DeterminantF;
        
    //4.-Left Cauchy-Green tensor b
    Matrix LeftCauchyGreen = prod(F0,trans(F0));

    //6.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	// e= 0.5*(1-invbT*invb)   
	this->CalculateAlmansiStrain(LeftCauchyGreen,StrainVector);
      }
 
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( 3 );

 		
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
    double& DeterminantF0               = rValues.GetDeterminantF0();
    const double& DeterminantF          = rValues.GetDeterminantF();

    double detF0 = DeterminantF0 * DeterminantF;

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
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_Cauchy);  //Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    DeterminantF0         *= DeterminantF;
  }

  //************************************************************************************
  //************************************************************************************


  void HyperElastic3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
  {
	    
    this->CalculateMaterialResponsePK1 (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step  
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK1,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    DeterminantF0         *= DeterminantF;

  }

  //************************************************************************************
  //************************************************************************************

  
  void HyperElastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
  {		
   
    this->CalculateMaterialResponseKirchhoff (rValues);
  
    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step  
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_Kirchhoff,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    DeterminantF0         *= DeterminantF;
  }


  //************************************************************************************
  //************************************************************************************

  void HyperElastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
  {

   
    this->CalculateMaterialResponseCauchy (rValues);
  
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    DeterminantF0         *= DeterminantF;


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
    Matrix StressMatrix( 3, 3 );

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

    rStressVector = MathUtils<double>::StressTensorToVector(StressMatrix,rStressVector.size());
    
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
		
    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rdetF0,rLameLambda,rLameMu,
								IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }

	  	
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

    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rinvF,rdetF0,rLameLambda,rLameMu,
								IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }
    	  
  }

  
  //***********************CONSTITUTIVE MATRIX COMPONENTS*******************************
  //************************************************************************************


  double& HyperElastic3DLaw::ConstitutiveComponent(double & rCabcd,
						   const Matrix &rMatrixIC, 
						   const double &rdetF0, 
						   const double &rLameLambda, 
						   const double &rLameMu, 
						   const unsigned int& a, const unsigned int& b, 
						   const unsigned int& c, const unsigned int& d)
  {
	  
    //(J²-1)/2
    //double auxiliar1 =  rdetF0*rdetF0;
    //double auxiliar2 =  (rdetF0*rdetF0-1);

    //(ln(J))
    double auxiliar1 =  1.0/rdetF0;
    double auxiliar2 =  (2.0*std::log(rdetF0));


    //1.Elastic constitutive tensor component:
    rCabcd=(rLameLambda*auxiliar1*rMatrixIC(a,b)*rMatrixIC(c,d));
    rCabcd+=((2*rLameMu-rLameLambda*auxiliar2)*0.5*(rMatrixIC(a,c)*rMatrixIC(b,d)+rMatrixIC(a,d)*rMatrixIC(b,c)));

    return rCabcd;
  }



  //***********************CONSTITUTIVE MATRIX COMPONENTS PULL-BACK*********************
  //************************************************************************************


  double& HyperElastic3DLaw::ConstitutiveComponent(double & rCabcd,
						   const Matrix &rMatrixIC, 
						   const Matrix &rinvF,
						   const double &rdetF0, 
						   const double &rLameLambda, 
						   const double &rLameMu, 
						   const unsigned int& a, const unsigned int& b, 
						   const unsigned int& c, const unsigned int& d)

  {

    rCabcd = 0;
    double Cijkl=0;
   
    unsigned int dimension = 3;

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
      {	  
	for(unsigned int l=0; l<dimension; l++)
	  {	   	    
	    for(unsigned int k=0; k<dimension; k++)
	      {
		for(unsigned int i=0; i<dimension; i++)
		  {
		    //Cijkl
		    rCabcd +=rinvF(a,i)*rinvF(b,j)*rinvF(c,k)*rinvF(d,l)*ConstitutiveComponent(Cijkl,rMatrixIC,rdetF0,rLameLambda,rLameMu,i,j,k,l);
		  }
	      }
	  }
      }

    return rCabcd;

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
