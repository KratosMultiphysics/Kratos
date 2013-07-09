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
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LinearElastic2DLaw::LinearElastic2DLaw(const LinearElastic2DLaw& rOther)
  : ConstitutiveLaw()
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
    return( rValue ); 
  }

  Vector& LinearElastic2DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
  {
    return( rValue );
  }

  Matrix& LinearElastic2DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
  {
    return( rValue );
  }


  //***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************


  void LinearElastic2DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

  }

  void LinearElastic2DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

  }

  void LinearElastic2DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
				    const ProcessInfo& rCurrentProcessInfo )
  {

  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::InitializeMaterial( const Properties& props,
					      const GeometryType& geom,
					      const Vector& ShapeFunctionsValues )
  {
 
  }

  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::InitializeSolutionStep( const Properties& props,
						  const GeometryType& geom, //this is just to give the array of nodes
						  const Vector& ShapeFunctionsValues,
						  const ProcessInfo& CurrentProcessInfo)
  {

  }
		
  //************************************************************************************
  //************************************************************************************

	
  void LinearElastic2DLaw::FinalizeSolutionStep( const Properties& props,
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


  void  LinearElastic2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
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

    //0.- Voigt size
    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);

    //3.-Determinant of the Total Deformation Gradient
    detF0 *= detF;
   
    //4.-Right Cauchy Green
    Matrix RightCauchyGreen = prod(trans(F0),F0);

    //5.-Inverse of the Right Cauchy-Green tensor C:
    double Trace_C=0;
    Matrix InverseRightCauchyGreen ( 2 , 2 );
    MathUtils<double>::InvertMatrix( RightCauchyGreen, InverseRightCauchyGreen, Trace_C);

    //6.-Green-Lagrange Strain:
    if(Options.Is( COMPUTE_STRAIN ))
      {
	//E= 0.5*(FT*F-1)
	StrainVector[0] = 0.5 * ( RightCauchyGreen( 0, 0 ) - 1.00 );
	StrainVector[1] = 0.5 * ( RightCauchyGreen( 1, 1 ) - 1.00 );
	StrainVector[2] = RightCauchyGreen( 0, 1 );
      }

    //7.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( 2 );

   
    if( Options.Is( COMPUTE_STRESS ) ){
	  
      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

      CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );		

    }
    else if(  Options.IsNot( COMPUTE_STRESS ) && Options.Is( COMPUTE_CONSTITUTIVE_TENSOR ) ){

      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

    }

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Strain "<<StrainVector<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;
		
  }


  //************************************************************************************
  //************************************************************************************


  void LinearElastic2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
  {
    CalculateMaterialResponsePK2 (rValues);
  
    Vector& StressVector=rValues.GetStressVector();
    const Matrix& F     =rValues.GetDeformationGradientF();
    const double& detF  =rValues.GetDeterminantF();
    
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
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   detF                  = rValues.GetDeterminantF(); 

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& detF0                         = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Voigt size
    Flags &Options=rValues.GetOptions();

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    //2.-Total Deformation Gradient
    Matrix F0 = prod(DeformationGradientF,DeformationGradientF0);

    //3.-Determinant of the Total Deformation Gradient
    detF0 *= detF;
        
    //4.-Left Cauchy-Green tensor b
    Matrix LeftCauchyGreen = prod(F0,trans(F0));

    //6.-Almansi Strain:
    if(Options.Is( COMPUTE_STRAIN ))
      {
	// e= 0.5*(1-invbT*invb)   
	Matrix InverseLeftCauchyGreen ( 2 , 2 );
	double Trace_b = 0;
	MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, Trace_b);

	Vector StrainVector( 3 );
	StrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
	StrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
	StrainVector[2] = InverseLeftCauchyGreen( 0, 1 );
       }
 
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( 2 );

    //7.-Incremental form
   
    if( Options.Is( COMPUTE_STRESS ) ){
	  
      CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

      CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );		

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
    double& detF0                       = rValues.GetDeterminantF0();;

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


  void LinearElastic2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
  {
    CalculateMaterialResponsePK1 (rValues);
  
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

  
  void LinearElastic2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
  {		
        CalculateMaterialResponseKirchhoff (rValues);
  
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

  void LinearElastic2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
  {
     
    CalculateMaterialResponseCauchy (rValues);
  
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& detF0                          =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();

    //2.-Update Internal Variables
    DeformationGradientF0  = prod(DeformationGradientF,DeformationGradientF0);
    detF0                  = MathUtils<double>::Det(DeformationGradientF0);


  }

  

  //***********************COMPUTE TOTAL STRESS PK2*************************************
  //************************************************************************************


  void LinearElastic2DLaw::CalculateStress( const Vector & rStrainVector,
					   const Matrix & rConstitutiveMatrix,
					   Vector& rStressVector )
  {
      
    //1.-2nd Piola Kirchhoff StressVector increment
    rStressVector = prod(rConstitutiveMatrix,rStrainVector);


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
