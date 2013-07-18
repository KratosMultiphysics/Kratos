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
#include "custom_constitutive/hyperelastic_mixed_2D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticMixed2DLaw::HyperElasticMixed2DLaw()
  : HyperElasticMixed3DLaw()
  {

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticMixed2DLaw::HyperElasticMixed2DLaw(const HyperElasticMixed2DLaw& rOther)
  : HyperElasticMixed3DLaw()
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElasticMixed2DLaw::Clone() const
  {
    HyperElasticMixed2DLaw::Pointer p_clone(new HyperElasticMixed2DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticMixed2DLaw::~HyperElasticMixed2DLaw()
  {
  }

  

  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElasticMixed2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
							Vector& rStrainVector )
  {

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

  }



  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElasticMixed2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						  Vector& rStrainVector )
  {

    // e= 0.5*(1-invbT*invb)   
    Matrix InverseLeftCauchyGreen ( 2 , 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = InverseLeftCauchyGreen( 0, 1 );


  }


  //***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
  //************************************************************************************


  void HyperElasticMixed2DLaw::CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
								     const Vector & rIsoStressVector,
								     const double & rdetF0,
								     const double & rTrace,
								     const double & rLameLambda,
								     const double & rLameMu,
								     Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };

    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<3; i++)
      {
	for(unsigned int j=0; j<3; j++)
	  {
	    rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rdetF0,rTrace,rLameLambda,rLameMu,IsoStressMatrix,
									 IndexVoigt2D[i][0],IndexVoigt2D[i][1],IndexVoigt2D[j][0],IndexVoigt2D[j][1]);
	  }

      }

	  	
  }

  //***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
  //************************************************************************************


  void HyperElasticMixed2DLaw::CalculateVolumetricConstitutiveMatrix (const Matrix & rMatrixIC,
								      const double & rdetF0,
								      const double & rLameLambda,
								      const double & rLameMu,
								      const GeometryType& rDomainGeometry,
								      const Vector & rShapeFunctions,
								      Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };	
    
    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rDomainGeometry, rShapeFunctions, Pressure);

    for(unsigned int i=0; i<3; i++)
      {
	for(unsigned int j=0; j<3; j++)
	  {
	    rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rdetF0,rLameLambda,rLameMu,Pressure,
									  IndexVoigt2D[i][0],IndexVoigt2D[i][1],IndexVoigt2D[j][0],IndexVoigt2D[j][1]);
	  }

      }

	  	
  }


  //******************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX PUSH-FORWARD****************
  //************************************************************************************

  void HyperElasticMixed2DLaw::CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
								     const Vector & rIsoStressVector,
								     const Matrix & rF,
								     const double & rdetF0,
								     const double & rTrace,
								     const double & rLameLambda,
								     const double & rLameMu,
								     Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };		
    
    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<3; i++)
      {
	for(unsigned int j=0; j<3; j++)
	  {
	    rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rF,rdetF0,rTrace,rLameLambda,rLameMu,IsoStressMatrix,
									 IndexVoigt2D[i][0],IndexVoigt2D[i][1],IndexVoigt2D[j][0],IndexVoigt2D[j][1]);
	  }

      }

	  	
  }

  //****************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX PUSH-FORWARD*****************
  //************************************************************************************

  void HyperElasticMixed2DLaw::CalculateVolumetricConstitutiveMatrix (const Matrix & rMatrixIC,
								      const Matrix & rF,
								      const double & rdetF0,
								      const double & rLameLambda,
								      const double & rLameMu,
								      const GeometryType& rDomainGeometry,
								      const Vector & rShapeFunctions,
								      Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int IndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };		
    
    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rDomainGeometry, rShapeFunctions, Pressure);


    for(unsigned int i=0; i<3; i++)
      {
    	for(unsigned int j=0; j<3; j++)
    	  {
	    rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ),rMatrixIC,rF,rdetF0,rLameLambda,rLameMu,Pressure,
									  IndexVoigt2D[i][0],IndexVoigt2D[i][1],IndexVoigt2D[j][0],IndexVoigt2D[j][1]);
	  }

      }

	  	
  }



} // Namespace Kratos
