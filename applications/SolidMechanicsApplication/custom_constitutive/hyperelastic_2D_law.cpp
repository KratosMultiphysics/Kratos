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
  : HyperElastic3DLaw()
  {

  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElastic2DLaw::HyperElastic2DLaw(const HyperElastic2DLaw& rOther)
  : HyperElastic3DLaw()
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

  

  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
							Vector& rStrainVector )
  {

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );



  }



  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rConstitutiveMatrix)
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


  //**************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX PULL-BACK*********************
  //************************************************************************************

  void HyperElastic2DLaw::CalculateConstitutiveMatrix (const Matrix &rMatrixIC,
						       const Matrix &rinvF,
						       const double &rdetF0,
						       const double &rLameLambda,
						       const double &rLameMu,
						       Matrix& rConstitutiveMatrix)
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



} // Namespace Kratos
