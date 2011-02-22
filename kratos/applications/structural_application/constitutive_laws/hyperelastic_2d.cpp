
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
*   Date:                $Date: 2010-02-05 16:10:12 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

/// BELOW ///////// TANGENT MATRIX

// System includes 
#include <iostream>
#include <limits>
// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/hyperelastic_2d.h" // USAR 2D com "D" maiusculo?

#include "includes/constitutive_law.h"  // nelson
#include "custom_utilities/tensor_utils.h"  // nelson
#include "custom_utilities/sd_math_utils.h"  // nelson

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{	/// HYPERELASTIC2D 
	/**Defines a hyperelastic constitutive law in 2D space, 
	for incompressible isotropic hyperelastic materials.
	//////////////////
	*/

    namespace Hyperelastic2DAuxiliaries
    {
//         boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
//         #pragma omp threadprivate(mstemp)
//         boost::numeric::ublas::bounded_matrix<double,2,2> msaux;
//         #pragma omp threadprivate(msaux)
	//	/// NEW ///
// 		Vector E1; // variables for pertubation method
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(E1)
// 		#endif
// 		Vector S1;
//                 //array_1d<double,6> S1;	
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(S1)
// 		#endif
// 		Vector E2; 
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(E2)
// 		#endif
// 		Vector S2;
//                 //array_1d<double,6> S2;
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(S2)
// 		#endif


    }
    using namespace Hyperelastic2DAuxiliaries;


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
		return false;
	}
	
	bool Hyperelastic2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
    void Hyperelastic2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Hyperelastic2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Hyperelastic2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void Hyperelastic2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
    }

      // NEW GET VALUE !!!!!!!!!
       double Hyperelastic2D::GetValue( const Variable<double>& rThisVariable )
  	{
  	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
  	}
	Vector& Hyperelastic2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
        if( rThisVariable == INTERNAL_VARIABLES )
        {
            rValue = ZeroVector(1);            
            return( rValue );
        }
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Hyperelastic2D::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
	}
    // NEW GET VALUE !!!!!!!!!

// 	bool Hyperelastic2D::Has( const Variable<double>& rThisVariable )
// 	{
// 		return false;
// 	}
// 	
// 	bool Hyperelastic2D::Has( const Variable<Vector>& rThisVariable )
// 	{
// 		if( rThisVariable == INSITU_STRESS ) // I think I need to change it to PK2_STRESS_TENSOR
// 			return true;
// 		return false;
// 	}
// 	
// 	bool Hyperelastic2D::Has( const Variable<Matrix>& rThisVariable )
// 	{
// 		return false;
// 	}
// 	
// 	double Hyperelastic2D::GetValue( const Variable<double>& rThisVariable )
// 	{
// 	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
// 	}
// 	
// 	Vector& Hyperelastic2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
// 	{
// 		if( rThisVariable == INSITU_STRESS )
// 			return mInSituStress;
// 		if( rThisVariable == INTERNAL_VARIABLES)
// 		  {
// 		    rValue=ZeroVector(1);
// 		    return (rValue);
// 		  }
// 	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
// 	}
// 	
// 	Matrix Hyperelastic2D::GetValue( const Variable<Matrix>& rThisVariable )
// 	{
// 	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
// 	}
// 
// 	void Hyperelastic2D::SetValue( const Variable<double>& rThisVariable, const double rValue, 
// 								   const ProcessInfo& rCurrentProcessInfo )
// 	{
// 	}
// 	
// 	void Hyperelastic2D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
// 								   const array_1d<double,3>& rValue, 
// 								   const ProcessInfo& rCurrentProcessInfo )
// 	{
// 	}
// 	
// 	void Hyperelastic2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
// 								 const ProcessInfo& rCurrentProcessInfo )
// 	{
// 		if( rThisVariable == INSITU_STRESS )
// 		{
// 			mInSituStress = rValue;
// 		}
// 	}
// 	
// 	void Hyperelastic2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
// 								 const ProcessInfo& rCurrentProcessInfo )
// 	{
// 	}
/*	
	void Hyperelastic2D::Calculate(const Variable<double>& rVariable, 
				      double& Output,
				      const ProcessInfo& rCurrentProcessInfo)
	  {
	    //Output = sqrt(mE/mDE);
	  }*/
	


	/**
	 *	TO BE TESTED!!!
	 */
	void Hyperelastic2D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
        //mCurrentStress = ZeroVector(6);
		mMU = props[MU]; // shear modulus
		mK = props[BULK_MODULUS]; // bulk modulus
 		mE = props[YOUNG_MODULUS];
 		mNU = props[POISSON_RATIO];
		mThickness = props[THICKNESS]; // to check
		mA0= geom.Area();  //geom[A0]; // to check
		mA= geom.Area();  // JUST A TEST, need to check if it will be updated later!!! 
		mAlpha=props[ALPHA];  // angle used for retraction
		mRetractionTime=props[RETRACTION_TIME];  // Time to start retraction
		/*mCtangent.resize(6,6,false);
		noalias(mCtangent) = ZeroMatrix(6,6);*/
		mMaterialParameters = props[MATERIAL_PARAMETERS];
		mInSituStress = ZeroVector(6);
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
    
    void Hyperelastic2D::ResetMaterial( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues)
    {
                  
    }

    void Hyperelastic2D::CalculateMaterialResponse( const Vector& StrainVector,
				      const Matrix& DeformationGradient,
				      Vector& StressVector, 
				      Matrix& AlgorithmicTangent,
				      const ProcessInfo& CurrentProcessInfo,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      bool CalculateStresses,
				      int CalculateTangent,
				      bool SaveInternalVariables)
    {
         	mA=geom.Area();//props[AREA]; // to check
 
		mcurrentThickness = (mA0*mThickness)/mA; // to check
		// RETRACTION ///////////////
		if (CurrentProcessInfo[TIME] >  mRetractionTime)
		//if (CurrentProcessInfo[TIME] > 0.1)   //0.01
		{
		mRetraction=mAlpha*CurrentProcessInfo[TIME]; // [DELTA_TIME]
		//KRATOS_WATCH("retraction");
		}
		else
		{
		mRetraction=0.0;
	        }                   

//CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);    
		CalculateStress(StrainVector, StressVector);
		CalculateConstitutiveMatrix(StrainVector, Ctang);  // nelson
		//CalculateConstitutiveMatrix(StrainVector, rResult);  // nelson
// 		mA=geom.Area();//props[AREA]; // to check
//  
// 		mcurrentThickness = (mA0*mThickness)/mA; // to check
// 		// RETRACTION ///////////////
// 		if (CurrentProcessInfo[TIME] > 0.0000001)   //0.01
// 		{
// 		mRetraction=mAlpha*CurrentProcessInfo[TIME]; // [DELTA_TIME]
// 		//KRATOS_WATCH("retraction");
// 		}
// 		else
// 		{
// 		mRetraction=0.0;
// 	        }
// 		//KRATOS_WATCH(CurrentProcessInfo[TIME]);
// 		//KRATOS_WATCH(mRetraction);
// 		// RETRACTION ///////////////
// 		//KRATOS_WATCH(mA0);
// 		//KRATOS_WATCH(mA);
// 		//KRATOS_WATCH(mThickness);
// 		//KRATOS_WATCH(mcurrentThickness);

    }
		
	/**
	 *	TO BE TESTED!!!
	 */
	
	// Implementing the Strain Energy for the Hyperelastic material (matrix)

	// this function gives W (strain energy), S (2nd Piola-Kirch) and CC (elasticity tensor) or Ctang (tangent matrix)
		// Ctang is also known as the Material Tangent Constitutive Tensor

	// ORDER:
		//#1 void Isotropic2D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
		//#2	void Isotropic2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
		//#3	void Isotropic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
		//#4	void Isotropic2D::CalculateCauchyStresses(Vector& rCauchy_StressVector,
		//		const Matrix& rF,
		//		const Vector& rPK2_StressVector,
		//		const Vector& rGreenLagrangeStrainVector)

	//void Hyperelastic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//*************************************************************************************************************************
	void Hyperelastic2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		// this funcion CalculateStress is called by EBST
		if( StressVector.size() != 3 ) // EBST
		{
			StressVector.resize(3,false); // false para no hacer copias
		}
		//noalias(StressVector) = prod(mCtangent,StrainVector);
		//mCurrentStress = StressVector;
		//noalias(StressVector) -= mInSituStress;


		////TRIAL
		//Hyperelastic2D::CalculateStressVector(StrainVector, auxStressVector);
		CalculateStressVector(StrainVector,auxStressVector);
		////TRIAL
		StressVector[0] = auxStressVector[0];
		StressVector[1] = auxStressVector[1];
		StressVector[2] = auxStressVector[3];
		//KRATOS_WATCH(StressVector)

	}
		// PS NOALIAS SOLO EN LADO ESQUIERDO!
//*************************************************************************************************************************
	/**
	Calculate Second Piola 
	*/

        void Hyperelastic2D::CalculateStressVector(const Vector& StrainVector, Vector& auxStressVector)
	//void Hyperelastic2D::CalculateStressVector(const Vector& StrainVector, array_1d<double,6>& rResult)
		{	
			// NEW // EBST
			if( auxStressVector.size() != 6 ) 
			{
			auxStressVector.resize(6,false); // false para no hacer copias
			}
			// NEW //

			Matrix I(3,3); // creating I matrix (Identity)
			//I=new IdentityMatrix(3);
			noalias(I)=ZeroMatrix(3,3);
			I(0,0)=1.0;
			I(1,1)=1.0;
			I(2,2)=1.0;
			
			Matrix C(3,3); // creating matrix C (C=Ft.F)
			noalias(C)=ZeroMatrix(3,3);
			
			// for EBST and Membrane
	      		if (mcurrentThickness == 0.0)  // DO I NEED THIS?!?!?! 
			{
			mcurrentThickness=(mA0/mA)*mThickness;
			}	
			//KRATOS_WATCH(mA0);
			//KRATOS_WATCH(mA);
			//KRATOS_WATCH(mThickness);
			//KRATOS_WATCH(mcurrentThickness);
			double Ez = 0.0;
			// Uncomment below if you want to use Ez
			Ez = (1.0/2.0)*((pow(mcurrentThickness,2.0)-pow(mThickness,2.0))/pow(mThickness,2.0)); 
			//KRATOS_WATCH(Ez);
			// TO INCLUDE RETRACTION !!!!!!!!!!
			Vector auxStrainVector(6); // 6 components
			auxStrainVector[0] = StrainVector[0]-mRetraction*(0.8*StrainVector[0]);
			auxStrainVector[1] = StrainVector[1]-mRetraction*(0.8*StrainVector[1]);
			auxStrainVector[2] = Ez;
			auxStrainVector[3] = StrainVector[2];
			auxStrainVector[4] = 0.0;
			auxStrainVector[5] = 0.0;
			//KRATOS_WATCH(auxStrainVector);
			//Vector auxStressVector(6); // 6 components
			
			//// IN TOTAL LAGRANGIAN! (BELOW)
			C(0,0) = 2.0 * auxStrainVector[0] + 1.0;
			C(1,1) = 2.0 * auxStrainVector[1] + 1.0;
			C(2,2) = 2.0 * auxStrainVector[2] + 1.0;
			C(0,1)=auxStrainVector[3];
			C(1,2)=auxStrainVector[4];
			C(0,2)=auxStrainVector[5]; 
			C(1,0)=C(0,1);
			C(2,1)=C(1,2);
			C(2,0)=C(0,2);
			//KRATOS_WATCH(C)
	
			Matrix Cinv(3,3); 
			noalias(Cinv)=ZeroMatrix(3,3);
			double det=0.0; // det of C  ( det(C)=J2 )
			//MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed

			double J2=MathUtils<double>::Det3( C ); // determinant of C is computed
			//KRATOS_WATCH(J2);
			double J=0.0;
			J=sqrt(J2); // J is computed
			//KRATOS_WATCH(J);


		// computing the derivative of strain energy: dW = S (the deviatoric part)
			
			double W=0.0;
			// Vector dW(3);
			Matrix dW(3,3);
			noalias(dW)=ZeroMatrix(3,3);
			double trC=0.0;
			trC=C(0,0)+C(1,1)+C(2,2); // computing Ckk (Trace of C = Sum of the diagonal terms) = I1 (First Invariant)
			//W=0.5*mMU*((pow(J2,0.3333333333333333333333333333))*trC-3.0);
			W=0.5*mMU*((1.0/(pow(J2,0.3333333333333333333333333333)))*trC-3.0); // ENERGY

			Matrix mdwtemp(3,3);
			noalias(mdwtemp)=ZeroMatrix(3,3);
			noalias(mdwtemp)=Cinv;
			mdwtemp*=trC;
			mdwtemp*=(-1.0/3.0);
			mdwtemp += I;	// parentesis solved  
			mdwtemp *= mMU;
			mdwtemp *= 1.0/pow(J2,0.3333333333333333333333333333);
			noalias(dW) = mdwtemp; // S computed!!! :)
			// Second Piola Kirchhoff is 3x3
			
			// puting second Piola Kirchhodd in Vector form
			Vector SPKir(6); // 6 components
			SD_MathUtils<double>::TensorToVector(dW,SPKir); // Rossi's Vector form
				
			//***************************

				
				Matrix Svol(3,3);
				noalias(Svol)=ZeroMatrix(3,3);
				double beta=9.0; // previous value: Beta=9
				double auxSvol=0.0;
				auxSvol=(mK/beta)*(1.0-(1.0/pow(J2,4.5))); // previous one! 
				//KRATOS_WATCH(Cinv);
				noalias(Svol)=Cinv;
				Svol*=auxSvol;
				//KRATOS_WATCH(Svol);
				Vector vSvol(6); // to put Svol in Vector form
				SD_MathUtils<double>::TensorToVector(Svol,vSvol); // Rossi's Vector form
				// Comment if not considering volumetric part
				//KRATOS_WATCH(SPKir);
				//KRATOS_WATCH(vSvol);
				noalias(auxStressVector)=SPKir+vSvol; 
				//KRATOS_WATCH(auxStressVector);
				//StressVector[0] = auxStressVector[0];
				//StressVector[1] = auxStressVector[1];
				//StressVector[2] = auxStressVector[3]; // xy component
				//KRATOS_WATCH(StressVector);	
						
			}
	/** 
	Calculate Constitutive Matrix, using the pertubation method
	*/
	void Hyperelastic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)

		{
			if( rResult.size1() != 3 ) // EBST
			{
			rResult.resize(3,3);
			}
// 
			 Vector E1(6);
			 Vector S1(6);	
			 Vector E2(6);
			 Vector S2(6);
			 // Vector E1;
// 			array_1d<double,6> S1;	
// 			Vector E2; 
// 			array_1d<double,6> S2;


		///////////// TANGENT MATRIX /////////// SERGIO OLLER + NELSON
	
		// PERTUBATION METHOD
			long double dE=0.00; // delta strain
			//double epsilon_real=std::numeric_limits<double>::epsilon(); //pertubation
			//KRATOS_WATCH(epsilon_real);
			long double epsilon=1E-9; // pertubation
			long double max=1E-14;
			E1.resize(6,false);
			//S1.resize(6,false);
			E2.resize(6,false);
			//S2.resize(6,false);

			// for EBST and Membrane
	      		mcurrentThickness=(mA0/mA)*mThickness;
			//KRATOS_WATCH(mA0);
			//KRATOS_WATCH(mA);
			//KRATOS_WATCH(mThickness);
			//KRATOS_WATCH(mcurrentThickness);
			double Ez = 0.0;
			// Uncomment below if you want to use Ez (uncomment the part related to y too! - line 587) 
			Ez = (1.0/2.0)*((pow(mcurrentThickness,2.0)-pow(mThickness,2.0))/pow(mThickness,2.0)); 
			//KRATOS_WATCH(Ez);
			
			Vector auxStrainVector(6); // 6 components
			auxStrainVector[0] = StrainVector[0]-mRetraction*(0.8*StrainVector[0]);  // 0.8*StrainVector[0] represents E0
                        auxStrainVector[1] = StrainVector[1]-mRetraction*(0.8*StrainVector[1]);
			auxStrainVector[2] = Ez;
			auxStrainVector[3] = StrainVector[2];
			auxStrainVector[4] = 0.0;
			auxStrainVector[5] = 0.0;
			//KRATOS_WATCH (auxStrainVector);
			Vector auxStressVector(6); // 6 components

			E1=auxStrainVector;
			//KRATOS_WATCH (E1);
			E2=auxStrainVector;
			//KRATOS_WATCH (E2);
			Matrix Ctot(6,6);  //
			noalias(Ctot)=ZeroMatrix(6,6);  // before Ctang
			//KRATOS_WATCH (Ctot);

			//// NEW FORM
			dE=(*std::min_element(auxStrainVector.begin(),auxStrainVector.end()))*epsilon;
			if (dE==0.00)
				{
					dE=epsilon;
				}
			if (dE<max)
				{
					dE=max;
				}
			for (unsigned int i=0;i<6;i++)
			{
				E2(i)+=dE; //pertubed Strain E2
				Hyperelastic2D::CalculateStressVector(E1,S1); // computing stress S1
				Hyperelastic2D::CalculateStressVector(E2,S2); // computing pertubed stress S2
				noalias(S2)=S2-S1;
				noalias(S2)=S2/(dE);
				for (unsigned int j=0;j<6;j++)
				{
					Ctot(j,i)=S2(j); // Tangent Matrix
					
				}
				
				//E2.resize(6,false);
				//S2.resize(6,false);
                                noalias(E2) = ZeroVector(6);
				noalias(S2) = ZeroVector(6);		
				//noalias(E2) = StrainVector; 
			}
			//noalias(rResult)=Ctang; // 2D
			//KRATOS_WATCH (Ctot);
			//noalias(rResult)=Ctang; // 2D
			
			// TEST !!! 
			Matrix auxCtang2(4,4);
			for (unsigned int i=0;i<4;i++)
			{
				for (unsigned int j=0;j<4;j++)
				{
					auxCtang2(i,j)=Ctot(i,j); // aux Tangent Matrix 4X4
					
				}
			}
			//KRATOS_WATCH (auxCtang2);
			// below putting zeros in last line and column
			Matrix auxCtang3(4,4);
			noalias(auxCtang3)=auxCtang2;
			for (unsigned int i=0;i<4;i++)
			{
				auxCtang3(2,i)=auxCtang2(3,i);
				auxCtang3(3,i)=auxCtang2(2,i);
				auxCtang3(i,2)=auxCtang2(i,3);
				auxCtang3(i,3)=auxCtang2(i,2);
			}
 			//auxCtang3(2,2)=auxCtang2(3,3);
 			//auxCtang3(2,3)=auxCtang2(3,2);
 			//auxCtang3(3,2)=auxCtang2(2,3);
 			//auxCtang3(3,3)=auxCtang2(2,2);
			//KRATOS_WATCH (auxCtang3(3,3));
			// end - Matrix 4X4 with last line and column with zeros (ref Sigma_Z)

			
			Matrix v(3,3);
			for (unsigned int i=0;i<3;i++)
			{
				for (unsigned int j=0;j<3;j++)
				{
				v(i,j)=auxCtang3(i,j);
				}
			}
			//KRATOS_WATCH (v);

			Matrix w(3,1);
			for (unsigned int i=0;i<3;i++)
			{
				w(i,0)=auxCtang3(i,3);
			}
			//KRATOS_WATCH (w);

			Matrix x(1,3);
			for (unsigned int i=0;i<3;i++)
			{
				x(0,i)=auxCtang3(3,i);
			}
			//KRATOS_WATCH (x);
			
			double y=auxCtang3(3,3);
			//KRATOS_WATCH (y);

			Matrix auxCtang(3,3);
			noalias(auxCtang) += v;
			//KRATOS_WATCH (auxCtang);
			Matrix auxtest(3,3);
			noalias(auxtest) = prod(w,x);
			//KRATOS_WATCH (auxtest);
// 			if (y != 0.0)
// 			{
// 			auxtest *= (1.0)/y;
// 			}
// 			else 
// 			{
// 			y = 0.000000000000000000001;
// 			auxtest *= (1.0)/y;
// 			}
			// Uncomment below if considering Ez 
			auxtest *= (1.0)/y;
			//KRATOS_WATCH (auxtest);			
			noalias(auxCtang) -= auxtest;
			//KRATOS_WATCH (auxCtang);

			noalias(rResult)=auxCtang;
			//KRATOS_WATCH (rResult);
	}

	void Hyperelastic2D::CalculateStressAndTangentMatrix(Vector& StressVector,  // check what is this function for!
			const Vector& StrainVector,
			Matrix& algorithmicTangent)
	{
	   // CalculateConsitutiveMatrix(StrainVector, algorithmicTangent);
	}

	
	//*************************************************************************************************************************
	/**
	Calculate Cauchy Stresses
	*/
	void Hyperelastic2D::CalculateCauchyStresses(	
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector, // in my case I already have S in a matrix form! How to read it?
		const Vector& rGreenLagrangeStrainVector)
	{
		//KRATOS_WATCH (rPK2_StressVector);
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
		//KRATOS_WATCH("Second Piola en HYP2D")
		//KRATOS_WATCH (S);
		
		double J = MathUtils<double>::Det2( rF );
		//KRATOS_WATCH("MATRIZ F em HYP2D")
		//KRATOS_WATCH (J);
		boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
		boost::numeric::ublas::bounded_matrix<double,2,2> msaux;

		
		noalias(mstemp) = prod(rF,S);
		//KRATOS_WATCH (mstemp);
		noalias(msaux) = prod(mstemp,trans(rF));
		//KRATOS_WATCH (msaux);
		msaux /= J;


		if(rCauchy_StressVector.size() != 3)
		rCauchy_StressVector.resize(3, false);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(0,1);
		//KRATOS_WATCH("CAUCHY IN HYPERELASTIC")
		//KRATOS_WATCH(rCauchy_StressVector)
	}

} // Namespace Kratos

//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************

/// UP ///////// TANGENT MATRIX


