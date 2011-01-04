/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
/* *********************************************************   
*          
*   Last Modified by:    $Author: virginia $
*   Date:                $Date: 2009-01-23 14:39:59 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/hyperelastic_2d.h" //OLD 3D

//#include "includes/constitutive_law.h"  DO I NEED TO INCLUDE THIS?!?!?!
//#include "C:\Kratos\kratosR1\applications\structural_application\custom_utilities/sd_math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{




	/**
	 *	TO BE TESTED!!!
	 */
	Hyperelastic2D::Hyperelastic2D() 
	: ConstitutiveLaw() 
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Hyperelastic2D::~Hyperelastic2D()
	{
	}
	
	bool Hyperelastic2D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Hyperelastic2D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS ) // I think I need to change it to PK2_STRESS_TENSOR
			return true;
		return false;
	}
	
	bool Hyperelastic2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& Hyperelastic2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
	    return( rValue );
	}
	
	Vector& Hyperelastic2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
	    return( rValue );
	}
	
	Matrix& Hyperelastic2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
	{
	    return( rValue );
	}

	void Hyperelastic2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hyperelastic2D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, //OLD 3
								   const array_1d<double,3>& rValue, //OLD 3
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hyperelastic2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void Hyperelastic2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Hyperelastic2D::InitializeMaterial( const Properties& props, //OLD 3D
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
        //mCurrentStress = ZeroVector(6);
		mMU = props[MU]; // Shear modulus
// 		mE = props[YOUNG_MODULUS];
// 		mNU = props[POISSON_RATIO];
		//KRATOS_WATCH(mMU);
		/*mCtangent.resize(6,6,false);
		noalias(mCtangent) = ZeroMatrix(6,6);*/
		
		mInSituStress = ZeroVector(3); //OLD 6
        //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);
//                 CalculateElasticMatrix(mCtangent, props[MATERIAL_PARAMETERS][0], props[MATERIAL_PARAMETERS][1]);
	}
	
	void Hyperelastic2D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void Hyperelastic2D::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			//mInSituStress -= mCurrentStress;
			//SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
		}
	}
    
    void Hyperelastic2D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {
		//KRATOS_TRY
        
		//CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);    
		CalculateConstitutiveMatrix(StrainVector, CC);  
		
		//KRATOS_CATCH(" ")
    }
		
	/**
	 *	TO BE TESTED!!!
	 */
	
	// Implementing the Strain Energy for the Hyperelastic material (matrix)

	// this function gives W (strain energy), S (2nd Piola-Kirch) and CC (elasticity tensor) 
	
	// OLD ONE!
	/*void Hyperelastic3D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)*/
	//

	// ORDER:
		//#1 void Isotropic3D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
		//#2	void Isotropic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
		//#3	void Isotropic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
		//#4	void Isotropic3D::CalculateCauchyStresses(Vector& rCauchy_StressVector,
		//		const Matrix& rF,
		//		const Vector& rPK2_StressVector,
		//		const Vector& rGreenLagrangeStrainVector)

	//void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//*************************************************************************************************************************
	void Hyperelastic2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 3 )  //OLD 6
		{
			StressVector.resize(3,false); // false para no hacer copias //OLD 6
		}
		//noalias(StressVector) = prod(mCtangent,StrainVector);
		//mCurrentStress = StressVector;
		//noalias(StressVector) -= mInSituStress;

		
			Matrix I(2,2); // creating I matrix (Identity) //OLD 3,3
			//I=new IdentityMatrix(3);
			noalias(I)=ZeroMatrix(2,2);//OLD 3,3
			I(0,0)=1.0;
			I(1,1)=1.0;
			//I(2,2)=1.0; //OLD : NOT COMMENTED
			
			Matrix C(2,2); // creating matrix C (C=Ft.F)//OLD 3,3
			noalias(C)=ZeroMatrix(2,2);//OLD 3,3

			//transforming "Strain Vector" into a symmetric 2nd order tensor "Voigt Notation"
			Matrix E(2,2);  //OLD 3,3
			/*Vector e(6);
			noalias(e)=StrainVector;*/
			//SD_MathUtils<double>::VectorToTensor(StrainVector,E);
			E=MathUtils<double>::StressVectorToTensor(StrainVector); // Rossi's Vector form
			//	MathUtils<double>::VectorToTensor(Vector& StrainVector,Matrix& E);

			 // transforming StrainVector into Matrix form
			Matrix auxE(2,2); //OLD 3,3
			noalias(auxE)=ZeroMatrix(2,2); //OLD 3,3
			//Matrix auxE=MathUtils<double>::StressVectorToTensor(StrainVector);
			noalias(auxE)=E;
			auxE *= 2;
			C += auxE;
			C += I; //NOW I HAVE "C" ( C = 2E+I )
				

			// The following is commented as I don't have rF as input:
			//double J = MathUtils<double>::Det3( rF ); // computing determinant of F
			//noalias(mstemp) = prod(trans(rF),rF); //computing C
		
			Matrix Cinv; 
			double det; // det of C  ( det(C)=J2 )
			//MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			MathUtils<double>::InvertMatrix2(C,Cinv,det); // Cinv and Det(C) are computed //OLD 3
			//MathUtils<double>::InvertMatrix(C,Cinv); 

			double J2=MathUtils<double>::Det2( C ); // determinant of C is computed //OLD 3
			double J;
			J=sqrt(J2); // J is computed



		// computing the derivative of strain energy: dW = S (the deviatoric part)
			
			double W;
			// Vector dW(3);
			Matrix dW(2,2); //OLD 3,3
			double trC;
			//trC=C(0,0)+C(1,1)+C(2,2); // computing Ckk (Trace of C = Sum of the diagonal terms) = I1 (First Invariant) //OLD
			trC=C(0,0)+C(1,1); // computing Ckk (Trace of C = Sum of the diagonal terms) = I1 (First Invariant) //NEW
			//W=0.5*mMU*((pow(J2,0.3333333333333333333333333333))*trC-3.0); OLD NOT COMMENTED
			
			// NEW BELOW
			double kG;
			kG=10.0;
			double G;
			G=(0.25)*(J2-1.0-2.0*log(J));
			// NEW: for 2D below
			W=0.5*mMU*trC-mMU*log(J)+kG*G;
			Vector S(3);
			Vector S1(3);
			Vector S2(3);
			Vector S3(3);
			Vector auxS1(3);
			auxS1(0)=1;
			auxS1(1)=1;
			auxS1(2)=0;
			S1=0.5*mMU*auxS1;
			S2(0)=0;
			S2(1)=0;
			S2(2)=0;
			SD_MathUtils<double>::TensorToVector(Cinv,S3);
			S3*=(-mMU/(2.0*J2)+kG/4.0-kG/4.0*J2)*J2;
			S=2.0*(S1+S2+S3);
			noalias(StressVector)=S;
			
			// END NEW

			//// OLD COMMENTED BELOW
			//Matrix mdwtemp(2,2);//OLD 3,3
			//mdwtemp=(-1/3)*trC*Cinv; 
			////Matrix I(3,3); // creating I matrix (Identity)
			////noalias(I)=ZeroMatrix(3,3);
			////I(0,0)=1.0;
			////I(1,1)=1.0;
			////I(2,2)=1.0;
			//mdwtemp += I;	// parentesis solved  
			//mdwtemp *= mMU;
			////mdwtemp *= (J^(2/3));
			//mdwtemp *= pow(J2,0.3333333333333333333333333333);
			//noalias(dW) = mdwtemp; // S computed!!! :)
			//// Second Piola Kirchhoff is 3x3
			//
			//// puting second Piola Kirchhodd in Vector form
			////Vector SPKir(6); // 6 components //OLD
			//Vector SPKir(3); // 3 components //NEW
			//SD_MathUtils<double>::TensorToVector(dW,SPKir); // Rossi's Vector form
			//// Uncomment below if not considering volumetric part
			////noalias(StressVector)=SPKir;
			////noalias(rPK2_StressVector)=noalias(SPKir); // do I need this step?
			//
			////***************************

			////// treating incompressibility !!!!!!!!!!!!!!!!!!!!
			//
			//	// Introducing Ogden model pg 244 - making it quasi-incompressible
			//	// volumetric Strain energy

			//	
			//	Matrix Svol(3,3);
			//	double k=1000000.0;
			//	double beta=9.0;
			//	double auxSvol;
			//	auxSvol=(k/beta)+(1-(1/pow(J,beta)));
			//	Svol=auxSvol*Cinv;
			//	//Vector vSvol(6); // to put Svol in Vector form //OLD
			//	Vector vSvol(3); // to put Svol in Vector form //NEW
			//	SD_MathUtils<double>::TensorToVector(Svol,vSvol); // Rossi's Vector form
			//	// Comment if not considering volumetric part
			//	noalias(StressVector)=SPKir+vSvol;
			//
			////// computing incompressibility !!!!!!!!!!!!!!!!!!!!

			//// using mixed elements:

			////	Vector EffectiveStress(6); // 6 components
			////	//EffectiveStress=J*Cinv;
			////	EffectiveStress*= J;
			////	EffectiveStress*= -1;
			////	EffectiveStress*= p; // need to define somewhere the value of p (Lagrange muliplier of the incompressibility equality constrain)
			////	//
			////***************************

			////mCurrentStress=StressVector; // What is this for?!?!?!? mCurrentStress???
			//
			//// Need to transform S in Sigma! How?!?! I need F !!! What the hell!!	
			////noalias(StressVector)=noalias(dW);
			// END OLD
	}
		// PS NOALIAS SOLO EN LADO ESQUIERDO!


//*************************************************************************************************************************
	void Hyperelastic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)

		{
			if( rResult.size1() != 3 )//OLD 6
			{
			rResult.resize(3,3,false);//OLD 6,6
			}
	
	// auxiliar computations

			Matrix I(2,2); // creating I matrix (Identity) //OLD 3,3
			//I=new IdentityMatrix(3);
			noalias(I)=ZeroMatrix(2,2);//OLD 3,3
			I(0,0)=1.0;
			I(1,1)=1.0;
			//I(2,2)=1.0;//OLD NOT COMMENTED
			
			Matrix C(2,2); // creating matrix C (C=Ft.F) //OLD 3,3
			noalias(C)=ZeroMatrix(2,2);//OLD 3,3

			//transforming "Strain Vector" into a symmetric 2nd order tensor "Voigt Notation"
			Matrix E(2,2);//OLD 3,3
			//SD_MathUtils<double>::VectorToTensor(StrainVector, E);// uncomment when using regular Voigt Notation
			E=MathUtils<double>::StressVectorToTensor(StrainVector); // using Rossi notation
			//MathUtils<double>::VectorToTensor(Vector& StrainVector,Matrix& E);

			 // transforming StrainVector into Matrix form
			Matrix auxE(2,2);//OLD 3,3
			noalias(auxE)=ZeroMatrix(2,2);//OLD 3,3
			//Matrix auxE=MathUtils<double>::StressVectorToTensor(StrainVector);
			noalias(auxE)=E;
			auxE *= 2;
			C += auxE;
			C += I; //NOW I HAVE "C" ( C = 2E+I )
			
			// The following is commented as I don't have rF as input:
			//double J = MathUtils<double>::Det3( rF ); // computing determinant of F
			//noalias(mstemp) = prod(trans(rF),rF); //computing C
		
			Matrix Cinv; 
			double det; // det of C  ( det(C)=J2 )
			//MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			MathUtils<double>::InvertMatrix2(C,Cinv,det); // Cinv and Det(C) are computed //OLD 3
			//MathUtils<double>::InvertMatrix(C,Cinv); 

			double J2=MathUtils<double>::Det2( C ); // determinant of C is computed //OLD 3

	// computing the Lagrangian (or material) elasticity tensor: CC

			Matrix CC(3,3);//OLD 6,6
			noalias(CC)=ZeroMatrix(3,3); // to check matrix size!//OLD 6,6
			
			//// Constant
			//double J23;
			////J23=1/(J^(2/3));
			//J23=1/pow(J2,0.3333333333333333333333333333); 
			//double constC; // 1st and 3rd terms. For the 2nd term, multiply by (-1/3)
			//constC=-(1/3)*mMU*J23;
			//// Volumetric part of CC
			//double constV2; // Constant related to the volumetric part of CC, 2nd term invC.invC
			//double constV3; // Constant related to the volumetric part of CC, 3rd term
			double kappa;
			kappa=10000.0; // bulk modulus (put a high value to simulate quasi-incompressibility)
			//constV2=(kappa/2)*((11/2)*(1/pow(J2,1.25))-(1/pow(J2,4.5)));
			//constV3=(kappa/2)*(1-(1/pow(J2,4.5)));


			// this matrix is used as the elastic matrix...to compute the Cauchy Stress!

				//	//filling the elasticity tensor CC
		

//////////////////////////////////////     START:       2D          /////////////////////////////////////////////////

			double const1;
			const1=kappa*J2;
			double const2;
			const2=kappa*(J2-1)-2*mMU;

			// Rossi notation 2D /// NEW !!!!!!!!!!!!!!!!!!!!!!

// 1ST ROW 
			CC(0,0) = (const1)*Cinv(0,0)*Cinv(0,0)+(const2)*(Cinv(0,0)*Cinv(0,0)+Cinv(0,0)*Cinv(0,0));
				CC(0,1) = (const1)*Cinv(0,0)*Cinv(1,1)+(const2)*(Cinv(0,1)*Cinv(0,1)+Cinv(0,1)*Cinv(0,1));  
					CC(0,2) = (const1)*Cinv(0,0)*Cinv(2,2)+(const2)*(Cinv(0,2)*Cinv(0,2)+Cinv(0,2)*Cinv(0,2));    

// 2ND ROW 
			CC(1,0) = CC(0,1);   
				CC(1,1) = (const1)*Cinv(1,1)*Cinv(1,1)+(const2)*(Cinv(1,1)*Cinv(1,1)+Cinv(1,1)*Cinv(1,1));
					 CC(1,2) = (const1)*Cinv(1,1)*Cinv(2,2)+(const2)*(Cinv(1,2)*Cinv(1,2)+Cinv(1,2)*Cinv(1,2)); 
					
// 3RD ROW 
			CC(2,0) = CC(0,2);		
				CC(2,1) = CC(1,2);    
					CC(2,2) = (const1)*Cinv(2,2)*Cinv(2,2)+(const2)*(Cinv(2,2)*Cinv(2,2)+Cinv(2,2)*Cinv(2,2));    


//////////////////////////////////////       END:     2D            /////////////////////////////////////////////////


	//// Rossi notation, with right constants

	//		// 1ST ROW  ----------> refference xx
	//		//CC(0,0) = constC*1constC(0)-((constC/3)+constV2)*2constC(0)+constC*;  
	//		 CC(0,0) = constC*I(0,0)*Cinv(0,0)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,0)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,0)+Cinv(0,0)*Cinv(0,0));
	//			CC(0,1) = constC*I(0,0)*Cinv(1,1)-((constC/3)+constV2)*Cinv(0,0)*Cinv(1,1)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(0,1)+Cinv(0,1)*Cinv(0,1));  
	//				CC(0,2) = constC*I(0,0)*Cinv(2,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(0,2)*Cinv(0,2)+Cinv(0,2)*Cinv(0,2));    
	//					CC(0,3) = constC*I(0,0)*Cinv(0,1)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,1)+Cinv(0,1)*Cinv(0,0));   
	//						CC(0,4) = constC*I(0,0)*Cinv(1,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(0,2)+Cinv(0,2)*Cinv(0,1));   
	//							CC(0,5) = constC*I(0,0)*Cinv(0,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,2)+Cinv(0,2)*Cinv(0,0));   
	//		// 2ND ROW  ----------> refference yy
	//			CC(1,0) = CC(0,1);   
	//				CC(1,1) = constC*I(1,1)*Cinv(1,1)-((constC/3)+constV2)*Cinv(1,1)*Cinv(1,1)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(1,1)+Cinv(1,1)*Cinv(1,1));
	//					 CC(1,2) = constC*I(1,1)*Cinv(2,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(1,2)*Cinv(1,2)+Cinv(1,2)*Cinv(1,2));   
	//						CC(1,3) = constC*I(1,1)*Cinv(0,1)-((constC/3)+constV2)*Cinv(1,1)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(1,0)*Cinv(1,1)+Cinv(1,1)*Cinv(1,0));   
	//							CC(1,4) = constC*I(1,1)*Cinv(1,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(1,2)+Cinv(1,2)*Cinv(1,1));   
	//								CC(1,5) = constC*I(1,1)*Cinv(0,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(1,0)*Cinv(1,2)+Cinv(1,2)*Cinv(1,0));   
	//		// 3RD ROW  ----------> refference zz
	//		CC(2,0) = CC(0,2);		
	//			CC(2,1) = CC(1,2);    
	//				CC(2,2) = constC*I(2,2)*Cinv(2,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(2,2)*Cinv(2,2)+Cinv(2,2)*Cinv(2,2));    
	//					CC(2,3) = constC*I(2,2)*Cinv(0,1)-((constC/3)+constV2)*Cinv(2,2)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(2,0)*Cinv(2,1)+Cinv(2,1)*Cinv(2,0)); 
	//						CC(2,4) = constC*I(2,2)*Cinv(1,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(2,1)*Cinv(2,2)+Cinv(2,1)*Cinv(2,2)); 
	//							CC(2,5) = constC*I(2,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(2,0)*Cinv(2,2)+Cinv(2,0)*Cinv(2,2)); 
	//		// 4TH ROW  ----------> refference to xy
	//		CC(3,0) = CC(0,3);   
	//			CC(3,1) = CC(1,3);   
	//				CC(3,2) = CC(2,3);      
	//					CC(3,3) = constC*I(0,1)*Cinv(0,1)-((constC/3)+constV2)*Cinv(0,1)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(1,1)+Cinv(0,1)*Cinv(1,0));    
	//						CC(3,4) = constC*I(0,1)*Cinv(1,2)-((constC/3)+constV2)*Cinv(0,1)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(1,2)+Cinv(0,2)*Cinv(1,1));
	//							CC(3,5) = constC*I(0,1)*Cinv(0,2)-((constC/3)+constV2)*Cinv(0,1)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(1,2)+Cinv(0,2)*Cinv(1,0));
	//		// 5TH ROW  ----------> refference to xz
	//		CC(4,0) = CC(0,4);   
	//			CC(4,1) = CC(1,4);   
	//				CC(4,2) = CC(2,4);   
	//					CC(4,3) = CC(3,4);   
	//						CC(4,4) = constC*I(1,2)*Cinv(1,2)-((constC/3)+constV2)*Cinv(1,2)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(2,1)); 
	//							CC(4,5) = constC*I(1,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(1,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(1,0)*Cinv(2,2)+Cinv(1,2)*Cinv(2,0)); 
	//		// 6TH ROW  ----------> refference to yz
	//		CC(5,0) = CC(0,5);   
	//			CC(5,1) = CC(1,5);   
	//				CC(5,2) = CC(2,5);   
	//					CC(5,3) = CC(3,5);   
	//						CC(5,4) = CC(4,5);    
	//							CC(5,5) = constC*I(0,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(0,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(2,2)+Cinv(0,2)*Cinv(2,0)); 
 //

// BELOW IS WORKING! 
				// Voigt Notation
						//					// 1ST ROW  ----------> refference xx
						////CC(0,0) = constC*1constC(0)-((constC/3)+constV2)*2constC(0)+constC*;  
						// CC(0,0) = constC*I(0,0)*Cinv(0,0)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,0)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,0)+Cinv(0,0)*Cinv(0,0));
						//	CC(0,1) = constC*I(0,0)*Cinv(1,1)-((constC/3)+constV2)*Cinv(0,0)*Cinv(1,1)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(0,1)+Cinv(0,1)*Cinv(0,1));  
						//		CC(0,2) = constC*I(0,0)*Cinv(2,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(0,2)*Cinv(0,2)+Cinv(0,2)*Cinv(0,2));    
						//			CC(0,3) = constC*I(0,0)*Cinv(1,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(0,2)+Cinv(0,1)*Cinv(0,2));   
						//				CC(0,4) = constC*I(0,0)*Cinv(0,2)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,2)+Cinv(0,0)*Cinv(0,2));   
						//					CC(0,5) = constC*I(0,0)*Cinv(0,1)-((constC/3)+constV2)*Cinv(0,0)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(0,1)+Cinv(0,0)*Cinv(0,1));   
						//// 2ND ROW  ----------> refference yy
						//	CC(1,0) = CC(0,1);   
						//		CC(1,1) = constC*I(1,1)*Cinv(1,1)-((constC/3)+constV2)*Cinv(1,1)*Cinv(1,1)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(1,1)+Cinv(1,1)*Cinv(1,1));
						//			 CC(1,2) = constC*I(1,1)*Cinv(2,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(1,2)*Cinv(1,2)+Cinv(1,2)*Cinv(1,2));   
						//				CC(1,3) = constC*I(1,1)*Cinv(1,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(1,2)+Cinv(1,1)*Cinv(1,2));   
						//					CC(1,4) = constC*I(1,1)*Cinv(0,2)-((constC/3)+constV2)*Cinv(1,1)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(1,0)*Cinv(1,2)+Cinv(1,0)*Cinv(1,2));   
						//						CC(1,5) = constC*I(1,1)*Cinv(0,1)-((constC/3)+constV2)*Cinv(1,1)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(1,0)*Cinv(1,1)+Cinv(1,0)*Cinv(1,1));   
						//// 3RD ROW  ----------> refference zz
						//CC(2,0) = CC(0,2);		
						//	CC(2,1) = CC(1,2);    
						//		CC(2,2) = constC*I(2,2)*Cinv(2,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(2,2)-((constC/2)+constV3)*(Cinv(2,2)*Cinv(2,2)+Cinv(2,2)*Cinv(2,2));    
						//			CC(2,3) = constC*I(2,2)*Cinv(1,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(2,1)*Cinv(2,2)+Cinv(2,1)*Cinv(2,2)); 
						//				CC(2,4) = constC*I(2,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(2,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(2,0)*Cinv(2,2)+Cinv(2,0)*Cinv(2,2)); 
						//					CC(2,5) = constC*I(2,2)*Cinv(0,1)-((constC/3)+constV2)*Cinv(2,2)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(2,0)*Cinv(2,1)+Cinv(2,0)*Cinv(2,1)); 
						//// 4TH ROW  ----------> refference to xy
						//CC(3,0) = CC(0,3);   
						//	CC(3,1) = CC(1,3);   
						//		CC(3,2) = CC(2,3);      
						//			CC(3,3) = constC*I(1,2)*Cinv(1,2)-((constC/3)+constV2)*Cinv(1,2)*Cinv(1,2)-((constC/2)+constV3)*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(1,2));    
						//				CC(3,4) = constC*I(1,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(1,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(2,2)+Cinv(1,2)*Cinv(0,2));
						//					CC(3,5) = constC*I(1,2)*Cinv(0,1)-((constC/3)+constV2)*Cinv(1,2)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,1)*Cinv(1,2)+Cinv(1,1)*Cinv(0,2));
						//// 5TH ROW  ----------> refference to xz
						//CC(4,0) = CC(0,4);   
						//	CC(4,1) = CC(1,4);   
						//		CC(4,2) = CC(2,4);   
						//			CC(4,3) = CC(3,4);   
						//				CC(4,4) = constC*I(0,2)*Cinv(0,2)-((constC/3)+constV2)*Cinv(0,2)*Cinv(0,2)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(2,2)+Cinv(0,2)*Cinv(0,2)); 
						//					CC(4,5) = constC*I(0,2)*Cinv(0,1)-((constC/3)+constV2)*Cinv(0,2)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(1,2)+Cinv(0,1)*Cinv(0,2)); 
						//// 6TH ROW  ----------> refference to yz
						//CC(5,0) = CC(0,5);   
						//	CC(5,1) = CC(1,5);   
						//		CC(5,2) = CC(2,5);   
						//			CC(5,3) = CC(3,5);   
						//				CC(5,4) = CC(4,5);    
						//					CC(5,5) = constC*I(0,1)*Cinv(0,1)-((constC/3)+constV2)*Cinv(0,1)*Cinv(0,1)-((constC/2)+constV3)*(Cinv(0,0)*Cinv(1,1)+Cinv(0,1)*Cinv(0,1)); 



		noalias(rResult)=CC;
	}
	//*************************************************************************************************************************
				  
	//void Hyperelastic3D::CalculateStress(const Vector& Strain Vector, Vector& Stress Vector)
	//{
	//	if (StressVector.size() !=6)
	//	{
	//		StressVector.resize(6);
	//	}
	//	noalias(StressVector)=prod(CC,StrainVector);
	//	mCurrentStress=StressVector;
	//	noalias(StressVector)-=mInSituStress;
	//}

	//void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	//{
	//	noalias(rResult)=noalias(CC);
	//}
	
	//*************************************************************************************************************************
	void Hyperelastic2D::CalculateCauchyStresses(	//OLD 3
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector, // in my case I already have S in a matrix form! How to read it?
		const Vector& rGreenLagrangeStrainVector)
	{
		//Matrix S=MathUtils<double>::StressVectorToTensor(rPK2_StressVector); // I think I don't need this!!! S is already a matrix!
		//Matrix S=dW;
		boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
		boost::numeric::ublas::bounded_matrix<double,2,2> msaux;
		
		Matrix S= MathUtils<double>::StressVectorToTensor(rPK2_StressVector);

		double J = MathUtils<double>::Det2( rF ); // computing determinant of F //OLD 3
		noalias(mstemp) = prod(rF,S); //computing temporary matrix F*S
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux *= J;

		if (rCauchy_StressVector.size()!=3) //OLD 6
		rCauchy_StressVector.resize(3);//OLD 6
		
		rCauchy_StressVector[0]=msaux(0,0);
		rCauchy_StressVector[1]=msaux(1,1);
		rCauchy_StressVector[2]=msaux(0,1); // NEW
		//rCauchy_StressVector[2]=msaux(2,2); // OLD NOT COMMENTED
		//rCauchy_StressVector[3]=msaux(1,2); // OLD NOT COMMENTED
		//rCauchy_StressVector[4]=msaux(1,3); // OLD NOT COMMENTED
		//rCauchy_StressVector[5]=msaux(2,3); // OLD NOT COMMENTED
	}

} // Namespace Kratos

	//*************************************************************************************************************************

	//////////////////////  !!!!!!!!!!!!! WHAT EVER !!!!!!!!!!!!!! //////////////////////////
//
//								/*noalias(mdwtemp)=prod(Cinv
//			noalias(mdwaux)=*/
//
//
//		// Matrix S = MathUtils<double>::StressVectorToTensor( dW );
//
//		// Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
//
//								////	Matrix::MatVec(EffectiveStress,mCtangent,StrainVector); 
//								//noalias(mEw) = prod(Cinv,EffectiveStress); 
//								////	Matrix::MatVec(mEw,Cinv,EffectiveStress); 
//								//	mEw -= StrainVector; 
//								//noalias(miaux) = prod(mstemp,trans(rF)); // where do I define "miaux"?
//								//
//								//noalias(msaux) = prod(mstemp,trans(rF));
//								//msaux *= J;
//
//								//if(rCauchy_StressVector.size() != 6)
//								//	rCauchy_StressVector.resize(6);
//								//
//								//rCauchy_StressVector[0] = msaux(0,0);
//								//rCauchy_StressVector[1] = msaux(1,1);
//								//rCauchy_StressVector[2] = msaux(2,2);
//								//rCauchy_StressVector[3] = msaux(1,2);
//								//rCauchy_StressVector[4] = msaux(1,3);
//								//rCauchy_StressVector[5] = msaux(2,3);
//	
//
//		//  void Hyperelastic3D::CalculateCauchyStresses(	Vector& rCauchy_StressVector, 
//		//												const Matrix& rF,	const Vector& rPK2_StressVector,	const Vector& rGreenLagrangeStrainVector)
//	
//	
//	//void Hyperelastic3D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
//	//{ 
//	//	//setting up material matrix
//	//	double c1 = E / ((1.00+NU)*(1-2*NU));
//	//	double c2 = c1 * (1-NU);
//	//	double c3 = c1 * NU;
//	//	double c4 = c1 * 0.5 * (1 - 2*NU);
//	//	//filling material matrix
//	//	C(0,0) = c2;    C(0,1) = c3;    C(0,2) = c3;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
//	//	C(1,0) = c3;    C(1,1) = c2;    C(1,2) = c3;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
//	//	C(2,0) = c3;    C(2,1) = c3;    C(2,2) = c2;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
//	//	C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = c4;    C(3,4) = 0.0;   C(3,5) = 0.0;
//	//	C(4,0) = 0.0;   C(4,1) = 0.0;   C(4,2) = 0.0;   C(4,3) = 0.0;   C(4,4) = c4;    C(4,5) = 0.0;
//	//	C(5,0) = 0.0;   C(5,1) = 0.0;   C(5,2) = 0.0;   C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = c4;
//	//	
//	//}
//	// TRY TO USE THE 2!!!! FOLLOWING STATEMENTS (Matrix and Vector forms of the Stress!)
//
//		 //virtual void CalculateStressAndTangentMatrix( Matrix& StressTensor,
//		 //                   const Matrix& StrainTensor,
//		 //                   MaterialTensorType& algorithmicTangent)
//		 //           {
//		 //               KRATOS_ERROR(std::logic_error,  
//		 //                            "Called the virtual function for CalculateStressAndTangentMatrix" , "");
//		 //           }
//		 //           
//		 //           /**
//		 //            * calculates the current stress and the material tangent
//		 //            * NOTE: there are two versions of this function: one for a matrix representation
//		 //            * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
//		 //  !!!!!!!!!!!!!!!!!!!!!!!!!          * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions 
//		 //            * available in MathUtils for either of them)
//		 //            * @param StressVector the calculated stress vector 
//		 //            * @param StrainVector the given strain vector
//		 //            * @param algorithmicTangent the calculated algorithmic tangent matrix
//		 //            */
//
//
//	//virtual void CalculateStressAndTangentMatrix( Vector& StressVector,
//		//                    const Vector& StrainVector,
//		//                    Matrix& algorithmicTangent)
//		//            {
//		//                KRATOS_ERROR(std::logic_error,  
//		//                             "Called the virtual function for CalculateStressAndTangentMatrix" , "");
//		//            }
//
//            /**
//             * this is to calculate the Cauchy stresses from a given Piola-Kirchhoff-2 stresses
//             * vector
//             * @param Cauchy_StressVector output: Cauchy stresses
//             * @param F input: deformation gradient
//             * @param PK2_StressVector input: Piola-Kirchhoff-2 stresses
//             * @param GreenLagrangeStrainVector input Green-Lagrange strains
//             */
//
//
//
//	/**
//	 *	TO BE TESTED!!!
//	 */
//	void Hyperelastic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
//	{
//		if( StressVector.size() != 6 )
//		{
//			StressVector.resize(6);
//		} 
//		noalias(StressVector) = prod(mCtangent,StrainVector);
//		mCurrentStress = StressVector;
//		noalias(StressVector) -= mInSituStress;
//// 		double c1 = mE / ((1.00+mNU)*(1-2*mNU));
//// 		double c2 = c1 * (1-mNU);
//// 		double c3 = c1 * mNU;
//// 		double c4 = c1 * 0.5 * (1 - 2*mNU);
//// 
//// 		StressVector[0] = c2*StrainVector[0] + c3 * (StrainVector[1] + StrainVector[2])	;
//// 		StressVector[1] = c2*StrainVector[1] + c3 * (StrainVector[0] + StrainVector[2])	;
//// 		StressVector[2] = c2*StrainVector[2] + c3 * (StrainVector[0] + StrainVector[1])	;
//// 		StressVector[3] = c4*StrainVector[3];
//// 		StressVector[4] = c4*StrainVector[4];
//// 		StressVector[5] = c4*StrainVector[5];
//	}
//	
//	/**
//	 *	TO BE REVIEWED!!!
//	 */
//	void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//	{
//		rResult = mCtangent;
//	}
//	
//    //**********************************************************************
//    void Hyperelastic3D::CalculateCauchyStresses(
//		Vector& rCauchy_StressVector,
//		const Matrix& rF,
//		const Vector& rPK2_StressVector,
//		const Vector& rGreenLagrangeStrainVector)
//    {
//		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
//
//        double J = MathUtils<double>::Det3( rF );
//
//		noalias(mstemp) = prod(rF,S);
//		noalias(msaux) = prod(mstemp,trans(rF));
//		msaux *= J;
//
//		if(rCauchy_StressVector.size() != 6)
//			rCauchy_StressVector.resize(6);
//		
//		rCauchy_StressVector[0] = msaux(0,0);
//		rCauchy_StressVector[1] = msaux(1,1);
//		rCauchy_StressVector[2] = msaux(2,2);
//		rCauchy_StressVector[3] = msaux(1,2);
//		rCauchy_StressVector[4] = msaux(1,3);
//		rCauchy_StressVector[5] = msaux(2,3);
//    }
//} // Namespace Kratos
