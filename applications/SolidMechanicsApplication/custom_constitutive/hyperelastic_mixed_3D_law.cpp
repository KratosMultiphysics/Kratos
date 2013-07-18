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
#include "custom_constitutive/hyperelastic_mixed_3D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticMixed3DLaw::HyperElasticMixed3DLaw()
  : HyperElastic3DLaw()
  {

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticMixed3DLaw::HyperElasticMixed3DLaw(const HyperElasticMixed3DLaw& rOther)
  : HyperElastic3DLaw()
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElasticMixed3DLaw::Clone() const
  {
    HyperElasticMixed3DLaw::Pointer p_clone(new HyperElasticMixed3DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticMixed3DLaw::~HyperElasticMixed3DLaw()
  {
  }

  

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*****************************MATERIAL RESPONSES*************************************
  //************************************************************************************


  void  HyperElasticMixed3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
  {

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF(); 

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Initialize parameters
    
    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

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

    SplitStressVector.Isochoric = ZeroVector(voigtsize);

    //OPTION 1: ( initial configuration )
    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( InverseRightCauchyGreen, IdentityMatrix, detF0, LameLambda, Trace_C, LameMu, StressMeasure_PK2, SplitStressVector.Isochoric );

    Vector IsochoricStressVector = SplitStressVector.Isochoric;

    
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

      SplitStressVector.Volumetric = ZeroVector(voigtsize);

      this->CalculateVolumetricStress ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

      //PK2 Stress:
      StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

      if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	StressVector = SplitStressVector.Isochoric;
      }
      else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	StressVector = SplitStressVector.Volumetric;
      }

    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
     
      this->CalculateIsochoricConstitutiveMatrix ( InverseRightCauchyGreen, IsochoricStressVector, detF0, Trace_C, LameLambda, LameMu, SplitConstitutiveMatrix.Isochoric );

      this->CalculateVolumetricConstitutiveMatrix ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

      //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
      ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric;
      
      if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric;
      }
      else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
      }
    
    }

    //OPTION 2: ( last known configuration : updated lagrangian approach only )
    if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) ){


      if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
	this->CalculateIsochoricStress( InverseRightCauchyGreen, IdentityMatrix, detF0, LameLambda, Trace_C, LameMu, StressMeasure_PK2, SplitStressVector.Isochoric );

      Vector IsochoricStressVector = SplitStressVector.Isochoric;
            
      TransformStresses(SplitStressVector.Isochoric,DeformationGradientF0,DeterminantF0,StressMeasure_PK2,StressMeasure_Cauchy);  //Cauchy Stress last known configuration == PK2 stress
      
      
      if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

	SplitStressVector.Volumetric = ZeroVector(voigtsize);

	this->CalculateVolumetricStress ( InverseRightCauchyGreen, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

	TransformStresses(SplitStressVector.Volumetric,DeformationGradientF0,DeterminantF0,StressMeasure_PK2,StressMeasure_Cauchy);  //Cauchy Stress last known configuration == PK2

	//PK2 Stress in the last known configuration:
	StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

	if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	  StressVector = SplitStressVector.Isochoric;
	}
	else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	  StressVector = SplitStressVector.Volumetric;
	}

      }

      if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
     
	this->CalculateIsochoricConstitutiveMatrix ( InverseRightCauchyGreen, IsochoricStressVector, DeformationGradientF0, detF0, Trace_C, LameLambda, LameMu, SplitConstitutiveMatrix.Isochoric );
	
	this->CalculateVolumetricConstitutiveMatrix ( InverseRightCauchyGreen, DeformationGradientF0, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

	SplitConstitutiveMatrix.Isochoric  /= DeterminantF0;
	SplitConstitutiveMatrix.Volumetric /= DeterminantF0;
	
	//if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
	ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric;
      
	if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	  ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric;
	}
	else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	  ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
	}
 
      }

    }
  
    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;
		
  }

  //************************************************************************************
  //************************************************************************************

  
  void HyperElasticMixed3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
  {

    //-----------------------------//
   
    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);
    
    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF(); 

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0(); 

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();      

    //-----------------------------//

    //0.- Initialize parameters
    
    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

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
        
    //4.-Left Cauchy-Green tensor b and Trace_b
    Matrix LeftCauchyGreen = prod(F0,trans(F0));

    double Trace_b = 0;
    for( unsigned int i=0; i<3; i++)
      {
	Trace_b += LeftCauchyGreen( i , i );
      }

    //6.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
      {
	// e= 0.5*(1-invbT*invb)   
	this->CalculateAlmansiStrain(LeftCauchyGreen,StrainVector);
      }
 
    //4.-Calculate Total PK2 stress   
    Matrix IdentityMatrix  = identity_matrix<double> ( 3 );
    
    SplitStressVector.Isochoric = ZeroVector(voigtsize);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( LeftCauchyGreen, IdentityMatrix, detF0, LameLambda, Trace_b, LameMu, StressMeasure_Kirchhoff, SplitStressVector.Isochoric );

    Vector IsochoricStressVector = SplitStressVector.Isochoric;

    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

      SplitStressVector.Volumetric = ZeroVector(voigtsize);

      this->CalculateVolumetricStress ( IdentityMatrix, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

      //Kirchhoff Stress:
      StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;
	

      if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	StressVector = SplitStressVector.Isochoric;
      }
      else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	StressVector = SplitStressVector.Volumetric;
      }

    }


    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){
     
      this->CalculateIsochoricConstitutiveMatrix ( IdentityMatrix, IsochoricStressVector, detF0, Trace_b, LameLambda, LameMu, SplitConstitutiveMatrix.Isochoric );

      this->CalculateVolumetricConstitutiveMatrix ( IdentityMatrix, detF0, LameLambda, LameMu, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

      //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
      ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric;
      
      if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ){    
	ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric;
      }
      else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ){
	ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
      }
    }

    
  }


  //******************************* COMPUTE ISOCHORIC STRESS  **************************
  //************************************************************************************



  void HyperElasticMixed3DLaw::CalculateIsochoricStress( const Matrix & rMatrixIC,
							 const Matrix & rIdentityMatrix,
							 const double & rdetF0,
							 const double & rTrace,
							 const double & rLameLambda, 
							 const double & rLameMu, 
							 StressMeasure rStressMeasure,
							 Vector& rIsoStressVector )
  {

    //1.-Identity build
    Matrix IsoStressMatrix ( 3, 3 );

    
    if(rStressMeasure == StressMeasure_PK2){    
      
      //rMatrixIC is InverseRightCauchyGreen

      //2.-Incompressible part of the 2nd Piola Kirchhoff Stress Matrix 
      IsoStressMatrix  = (rIdentityMatrix - (rTrace/3.0)*rMatrixIC );
      IsoStressMatrix *= rLameMu*pow(rdetF0,(-2.0/3.0));

      //std::cout<<" PK2 "<<std::endl;
    }
    
    if(rStressMeasure == StressMeasure_Kirchhoff){

      //rMatrixIC is LeftCauchyGreen
      
      //2.-Incompressible part of the Kirchhoff Stress Matrix 
      IsoStressMatrix  = (rMatrixIC - (rTrace/3.0)*rIdentityMatrix );
      IsoStressMatrix *= rLameMu*pow(rdetF0,(-2.0/3.0));

      //std::cout<<" Kirchooff "<<std::endl;

    }
 

    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());
    
  }

  //******************************* COMPUTE DOMAIN PRESSURE  ***************************
  //************************************************************************************


  double &  HyperElasticMixed3DLaw::CalculateDomainPressure (const GeometryType& rDomainGeometry,
							     const Vector & rShapeFunctions, 
							     double & rPressure)
  {

    const unsigned int number_of_nodes = rDomainGeometry.size();
    
    rPressure = 0;
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
 	rPressure += rShapeFunctions[j] * rDomainGeometry[j].GetSolutionStepValue(PRESSURE);
      }

    return rPressure;
  }


  //******************************* COMPUTE VOLUMETRIC STRESS  *************************
  //************************************************************************************

  void HyperElasticMixed3DLaw::CalculateVolumetricStress(const Matrix & rMatrixIC,
							 const double & rdetF0,
							 const double & rLameLambda, 
							 const double & rLameMu, 
							 const GeometryType& rDomainGeometry,
							 const Vector & rShapeFunctions,
							 Vector& rVolStressVector )
  {

    //1.- Declaration
    Matrix VolStressMatrix ( 3 , 3 );
    
    double Pressure = 0;

    Pressure = CalculateDomainPressure (rDomainGeometry, rShapeFunctions, Pressure);

    //2.- Volumetric part of the Kirchhoff StressMatrix from nodal pressures
    VolStressMatrix = rdetF0 * Pressure * rMatrixIC;
    

    //3.- Volumetric part of the  Kirchhoff StressMatrix

    /*
    //(J²-1)/2
    //double auxiliar1 =  0.5*(rdetF0*rdetF0-1);
    
    //(ln(J))
    double auxiliar1 =  std::log(rdetF0);
    
    double BulkModulus= rLameLambda + (2.0/3.0) * rLameMu;

    //2.-Volumetric part of the Kirchhoff Stress Matrix 

     VolStressMatrix  = BulkModulus * auxiliar1 * rMatrixIC; 
    */

    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());
    
  }

  //***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
  //************************************************************************************

  void HyperElasticMixed3DLaw::CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
								     const Vector & rIsoStressVector,
								     const double & rdetF0,
								     const double & rTrace,
								     const double & rLameLambda,
								     const double & rLameMu,
								     Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };


    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rdetF0,rTrace,rLameLambda,rLameMu,IsoStressMatrix,
									 IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }

	  	
  }

  //***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
  //************************************************************************************


  void HyperElasticMixed3DLaw::CalculateVolumetricConstitutiveMatrix (const Matrix & rMatrixIC,
								      const double & rdetF0,
								      const double & rLameLambda,
								      const double & rLameMu,
								      const GeometryType& rDomainGeometry,
								      const Vector & rShapeFunctions,
								      Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    
    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rDomainGeometry, rShapeFunctions, Pressure);


    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rdetF0,rLameLambda,rLameMu,Pressure,
									  IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }

	  	
  }


  //******************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX PUSH-FORWARD****************
  //************************************************************************************

  void HyperElasticMixed3DLaw::CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
								     const Vector & rIsoStressVector,
								     const Matrix & rF,
								     const double & rdetF0,
								     const double & rTrace,
								     const double & rLameLambda,
								     const double & rLameMu,
								     Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };


    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rF,rdetF0,rTrace,rLameLambda,rLameMu,IsoStressMatrix,
									 IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }

	  	
  }

  //***************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX PUSH-FORWARD******************
  //************************************************************************************

  void HyperElasticMixed3DLaw::CalculateVolumetricConstitutiveMatrix (const Matrix & rMatrixIC,
								      const Matrix & rF,
								      const double & rdetF0,
								      const double & rLameLambda,
								      const double & rLameMu,
								      const GeometryType& rDomainGeometry,
								      const Vector & rShapeFunctions,
								      Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    
    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rDomainGeometry, rShapeFunctions, Pressure);

    for(unsigned int i=0; i<6; i++)
      {
	for(unsigned int j=0; j<6; j++)
	  {
	    rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rF,rdetF0,rLameLambda,rLameMu,Pressure,
									  IndexVoigt3D[i][0],IndexVoigt3D[i][1],IndexVoigt3D[j][0],IndexVoigt3D[j][1]);
	  }

      }

	  	
  }

  
  //********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
  //************************************************************************************


  double& HyperElasticMixed3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
								const Matrix & rMatrixIC, 
								const double & rdetF0, 
								const double & rTrace, 
								const double & rLameLambda, 
								const double & rLameMu,
								const Matrix & rIsoStressMatrix,
								const unsigned int& a, const unsigned int& b, 
								const unsigned int& c, const unsigned int& d)
  {
	  
    //Isochoric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)    

    rCabcd  = (1.0/3.0)*(rMatrixIC(a,b)*rMatrixIC(c,d));
    rCabcd -= (0.5*(rMatrixIC(a,c)*rMatrixIC(b,d)+rMatrixIC(a,d)*rMatrixIC(b,c)));
    rCabcd *= pow(rdetF0,(-2.0/3.0))*rTrace*rLameMu;
    
    rCabcd += (rMatrixIC(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rMatrixIC(a,b));
    rCabcd *= (-2.0/3.0);

    return rCabcd;
  }


  //********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT************************
  //************************************************************************************


  double& HyperElasticMixed3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
								 const Matrix & rMatrixIC, 
								 const double & rdetF0, 
								 const double & rLameLambda, 
								 const double & rLameMu,
								 const double & rPressure,
								 const unsigned int& a, const unsigned int& b, 
								 const unsigned int& c, const unsigned int& d)
  {

    //Volumetric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)    

    rCabcd = (rMatrixIC(a,b)*rMatrixIC(c,d));
    rCabcd -= 2*(0.5*(rMatrixIC(a,c)*rMatrixIC(b,d)+rMatrixIC(a,d)*rMatrixIC(b,c)));

    rCabcd *= rPressure*rdetF0;


    /*
    double BulkModulus= rLameLambda + (2.0/3.0) * rLameMu;

    //(J²-1)/2
    //double auxiliar1 =  rdetF0*rdetF0;
    //double auxiliar2 =  (rdetF0*rdetF0-1);
    //double auxiliar3 =  BulkModulus;

    //(ln(J))
    double auxiliar1 =  1.0;
    double auxiliar2 =  (2.0*std::log(rdetF0));
    double auxiliar3 =  BulkModulus/rdetF0;

    //1.Volumetric Elastic constitutive tensor component:
    rCabcd = auxiliar1*(rMatrixIC(a,b)*rMatrixIC(c,d));
    rCabcd -= auxiliar2*(0.5*(rMatrixIC(a,c)*rMatrixIC(b,d)+rMatrixIC(a,d)*rMatrixIC(b,c)));
    rCabcd *= auxiliar3;
    */



    return rCabcd;
  }


  //**************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT PUSH-FORWARD******************
  //************************************************************************************


  double& HyperElasticMixed3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
								const Matrix & rMatrixIC, 
								const Matrix & rF,
								const double & rdetF0, 
								const double & rTrace, 
								const double & rLameLambda, 
								const double & rLameMu,
								const Matrix & rIsoStressMatrix,
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
		    rCabcd +=rF(a,i)*rF(b,j)*rF(c,k)*rF(d,l)*IsochoricConstitutiveComponent(Cijkl,rMatrixIC,rdetF0,rTrace,rLameLambda,rLameMu,rIsoStressMatrix,i,j,k,l);
		  }
	      }
	  }
      }

    return rCabcd;


  }


  //**************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT PUSH-FORWARD*****************
  //************************************************************************************

  double& HyperElasticMixed3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
								  const Matrix & rMatrixIC, 
								  const Matrix & rF,
								  const double & rdetF0, 
								  const double & rLameLambda, 
								  const double & rLameMu,
								  const double  & rPressure,
								  const unsigned int& a, const unsigned int& b, 
								  const unsigned int& c, const unsigned int& d)
  {



    rCabcd = 0;
    double Cijkl = 0;
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
		    rCabcd +=rF(a,i)*rF(b,j)*rF(c,k)*rF(d,l)*VolumetricConstitutiveComponent(Cijkl,rMatrixIC,rdetF0,rLameLambda,rLameMu,rPressure,i,j,k,l);
		  }
	      }
	  }
      }

    return rCabcd;

  }




    
} // Namespace Kratos
