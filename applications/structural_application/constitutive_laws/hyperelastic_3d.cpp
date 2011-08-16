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
#include "constitutive_laws/hyperelastic_3d.h"

#include "includes/constitutive_law.h"  // nelson
#include "custom_utilities/tensor_utils.h"  // nelson
#include "custom_utilities/sd_math_utils.h"  // nelson

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{	/// HYPERELASTIC3D 
	/** Defines a hyperelastic constitutive law in 3D space, 
	// for incompressible isotropic hyperelastic materials.
	*/

    namespace Hyperelastic3DAuxiliaries
    {
//         boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
//         #pragma omp threadprivate(mstemp)
//         boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
//         #pragma omp threadprivate(msaux)
// 		/// NEW ///
// 		Vector E1; // variables for pertubation method
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(E1)
// 		#endif
// 		Vector S1;
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(S1)
// 		#endif
// 		Vector E2; 
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(E2)
// 		#endif
// 		Vector S2;
// 		#ifdef _OPENMP
// 		#pragma omp threadprivate(S2)
// 		#endif


    }
    using namespace Hyperelastic3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	Hyperelastic3D::Hyperelastic3D() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Hyperelastic3D::~Hyperelastic3D()
	{
	}
	
	bool Hyperelastic3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Hyperelastic3D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS ) 
			return true;
		return false;

	}
	
	bool Hyperelastic3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
// 	double Hyperelastic3D::GetValue( const Variable<double>& rThisVariable )
// 	{
// 	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
// 	}
	
// 	Vector& Hyperelastic3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) 
// 	{
// 		if( rThisVariable == INSITU_STRESS )
// 			return mInSituStress;
// 		if (rThisVariable == INTERNAL_VARIABLES)
// 		{
// 		  rValue=ZeroVector(1);
// 		  return(rValue);
// 		}
// 	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
// 	}
	
// 	Matrix Hyperelastic3D::GetValue( const Variable<Matrix>& rThisVariable )
// 	{
// 	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
// 	}

	void Hyperelastic3D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hyperelastic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo  )
	{
	}
	
// 	void Hyperelastic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
// 								 const ProcessInfo& rCurrentProcessInfo )
// 	{
// 		if( rThisVariable == INSITU_STRESS )
// 		{
// 			mInSituStress = rValue;
// 		}
// 	}
	
	void Hyperelastic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	  // Output = sqrt(mE/mDE); // E is young modulus and DE is density
	}
	

	void Hyperelastic3D::Calculate( const Variable<double>& rVariable, double& Output,
								 const ProcessInfo& rCurrentProcessInfo )
	{ // use to consider internal variables
	}


	/**
	 *	TO BE TESTED!!!
	 */
	void Hyperelastic3D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
        //mCurrentStress = ZeroVector(6);
		//mMaterialParameters = props[MATERIAL_PARAMETERS];
		mMU = props[MU]; // shear modulus
		mK = props[BULK_MODULUS]; // bulk modulus
 		//mE = props[YOUNG_MODULUS];
 		//mNU = props[POISSON_RATIO];
		//mDE = props[DENSITY]; 

		/*mCtangent.resize(6,6,false);
		noalias(mCtangent) = ZeroMatrix(6,6);*/
		
		mInSituStress = ZeroVector(6);
        //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);
//                 CalculateElasticMatrix(mCtangent, props[MATERIAL_PARAMETERS][0], props[MATERIAL_PARAMETERS][1]);
	}
	
	void Hyperelastic3D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void Hyperelastic3D::FinalizeSolutionStep( const Properties& props,
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
    
    void Hyperelastic3D::ResetMaterial( const Properties& props,  // is the old UpdateMaterial
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues )
                                      
    {
        //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);    
	//	CalculateConstitutiveMatrix(StrainVector, Ctang);  // nelson
    }
	

    void Hyperelastic3D::CalculateMaterialResponse( const Vector& StrainVector,  // is the old UpdateMaterial
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
        //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);    
	CalculateStress(StrainVector, StressVector);	
	CalculateConstitutiveMatrix(StrainVector, Ctang);  // nelson
    }
	
	/**
	 *	TO BE TESTED!!!
	 */
	
	// Implementing the Strain Energy for the Hyperelastic material (matrix)

	// this function gives W (strain energy), S (2nd Piola-Kirch) and CC (elasticity tensor) or Ctang (tangent matrix)
		// Ctang is also known as the Material Tangent Constitutive Tensor

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
	void Hyperelastic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		// this funcion CalculateStress is called by total_lagrangian.
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6,false); // false para no hacer copias
		}
		//noalias(StressVector) = prod(mCtangent,StrainVector);
		//mCurrentStress = StressVector;
		//noalias(StressVector) -= mInSituStress;


		////TRIAL
		Hyperelastic3D::CalculateStressVector(StrainVector, StressVector);
		////TRIAL
		

	}
		// PS NOALIAS SOLO EN LADO ESQUIERDO!
//*************************************************************************************************************************
	/**
	Calculate Second Piola 
	*/

	void Hyperelastic3D::CalculateStressVector(const Vector& StrainVector, Vector& StressVector)
		{	
			// NEW //
			if( StressVector.size() != 6 )
			{
			StressVector.resize(6,false); // false para no hacer copias
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
			
			//// IN TOTAL LAGRANGIAN! (BELOW)
			C(0,0) = 2.0 * StrainVector[0] + 1.0;
			C(1,1) = 2.0 * StrainVector[1] + 1.0;
			C(2,2) = 2.0 * StrainVector[2] + 1.0;
			C(0,1)=StrainVector[3];
			C(1,2)=StrainVector[4];
			C(0,2)=StrainVector[5]; 
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
				noalias(StressVector)=SPKir+vSvol; 
			
			}


//*************************************************************************************************************************
////ISOTROPIC
//	void Hyperelastic3D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
//	{ 
//		//setting up material matrix
//		double c1 = E / ((1.00+NU)*(1-2*NU));
//		double c2 = c1 * (1-NU);
//		double c3 = c1 * NU;
//		double c4 = c1 * 0.5 * (1 - 2*NU);
//		//filling material matrix
//		C(0,0) = c2;    C(0,1) = c3;    C(0,2) = c3;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
//		C(1,0) = c3;    C(1,1) = c2;    C(1,2) = c3;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
//		C(2,0) = c3;    C(2,1) = c3;    C(2,2) = c2;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
//		C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = c4;    C(3,4) = 0.0;   C(3,5) = 0.0;
//		C(4,0) = 0.0;   C(4,1) = 0.0;   C(4,2) = 0.0;   C(4,3) = 0.0;   C(4,4) = c4;    C(4,5) = 0.0;
//		C(5,0) = 0.0;   C(5,1) = 0.0;   C(5,2) = 0.0;   C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = c4;
//		
//	}
//	void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//	{
//		rResult = mCtangent;
//	}
//	//ISOTROPIC

	/** 
	Calculate Constitutive Matrix, using the pertubation method
	*/

	void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)

		{
			if( rResult.size1() != 6 )
			{
			rResult.resize(6,6);
			}
	
	
			 Vector E1(6);
			 Vector S1(6);	
			 Vector E2(6);
			 Vector S2(6);

		///////////// TANGENT MATRIX /////////// SERGIO OLLER + NELSON
		////	Vector S1(6);
		////	Vector S2(6);
		////	if((c.size1()!=6) || c.size2() !=6)
		////	{
		////		c.resize(6,6);
		////	}
		////		if((Ctang.size1()!=6) || Ctang.size2() !=6)
		////	{
		////		Ctang.resize(6,6);
		////	}
		////	//noalias(c)=ZeroMatrix(6,6);
		////	CalculateStressVector (StrainVector, S1);
		////	//KRATOS_WATCH(S1);
		////	//private CalculateStressVector (const Vector& StrainVector, S1)
		////	int i;
		////		for (i=0; i<6; i++)
		////		{
		////			Vector Ep(StrainVector);
		////			double epsilon=std::numeric_limits<double>::epsilon();
		////			Ep[i]+=epsilon;
		////			//KRATOS_WATCH(epsilon);
		////			CalculateStressVector (Ep, S2);
		////			//KRATOS_WATCH(S2);
		////			int j;
		////			for  (j=0; j<6; j++)
		////			{				
		////				c(j,i)=(1.0/epsilon)*(S2[j]-S1[i]);
		////			}
		////			
		////		}
		////		noalias(Ctang)=c;
		//
		// PERTUBATION METHOD
			long double dE=0.00; // delta strain
			//double epsilon_real=std::numeric_limits<double>::epsilon(); //pertubation
			//KRATOS_WATCH(epsilon_real);
			long double epsilon=1E-9; // pertubation
			long double max=1E-14;
			E1.resize(6,false);
			S1.resize(6,false);
			E2.resize(6,false);
			S2.resize(6,false);

			E1=StrainVector;
			E2=StrainVector;

			//Matrix Ctang(6,6);  //
			noalias(Ctang)=ZeroMatrix(6,6);
			//KRATOS_WATCH (Ctang);

			//// NEW FORM
			dE=(*std::min_element(StrainVector.begin(),StrainVector.end()))*epsilon;
			if (dE==0.00)
				{
					dE=epsilon;
				}
			if (dE<max)
				{
					dE=max;
				}
			for (unsigned int i=0;i<E1.size();i++)
			{
				E2(i)+=dE; //pertubed Strain E2
				Hyperelastic3D::CalculateStressVector(E1,S1); // computing stress S1
				Hyperelastic3D::CalculateStressVector(E2,S2); // computing pertubed stress S2
				noalias(S2)=S2-S1;
				noalias(S2)=S2/(dE);
				for (unsigned int j=0;j<E1.size();j++)
				{
					Ctang(j,i)=S2(j); // Tangent Matrix
					
				}
				
				E2.resize(6,false);
				S2.resize(6,false);	
				noalias(E2)=StrainVector;
			}
			noalias(rResult)=Ctang;
			//KRATOS_WATCH (rResult);
	}

			// BELOW // OLD FORM C_TANG
			//for (unsigned int i=0;i<E1.size();i++)
			//{
			//	if (fabs(StrainVector(i))<1E-14)
			//	{
			//	dE=(*std::min_element(StrainVector.begin(),StrainVector.end()))*epsilon;
			//		if (dE==0.00)
			//		{
			//		dE=epsilon;
			//		}
			//	}
			//	else
			//	{
			//		dE=StrainVector(i)*epsilon;
			//	}
			//	 //i won't use this as my epsilon is already lower than max).
			//	if (dE<max)
			//	{
			//		dE=max;
			//	}
			//	
			//	E1(i)+=dE; //adding pertubation to E1
			//	E2(i)-=dE; //subtracting pertubation to E2
			//	//// ME
			//	//E1(i)+=0; //adding pertubation to E1
			//	//E2(i)+=dE; //subtracting pertubation to E2
			//	//// ME

			//	Hyperelastic3D::CalculateStressVector(E1,S1); // computing pertubed stress S1
			//	Hyperelastic3D::CalculateStressVector(E2,S2); // computing pertubed stress S2

			//	noalias(S1)=S1-S2; // computing delta stress
			//	noalias(S1)=S1/(2.00*dE); // average
			//	//noalias(S1)=S2-S1; // computing delta stress // ME
			//	//noalias(S1)=S1/(dE); // average

			//	for (unsigned int j=0;j<E1.size();j++)
			//	{
			//		Ctang(j,i)+=S1(j); // Tangent Matrix
			//		
			//	}
			//	//KRATOS_WATCH (i);
			//	//KRATOS_WATCH (Ctang);
			//	S1=StrainVector;
			//	S2=StrainVector;		
	//		
	//		
	//		
	//		}
	//	
	//noalias(rResult)=Ctang;
	////KRATOS_WATCH (rResult);

	//}
	
	void Hyperelastic3D::CalculateStressAndTangentMatrix(Vector& StressVector,  // check what is this function for!
			const Vector& StrainVector,
			Matrix& algorithmicTangent)
	{
	   // CalculateConsitutiveMatrix(StrainVector, algorithmicTangent);
	}



	//*************************************************************************************************************************
	/**
	Calculate Cauchy Stresses
	*/
	void Hyperelastic3D::CalculateCauchyStresses(	
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector, // in my case I already have S in a matrix form! How to read it?
		const Vector& rGreenLagrangeStrainVector)
	{
		//Matrix S=MathUtils<double>::StressVectorToTensor(rPK2_StressVector); // I think I don't need this!!! S is already a matrix!
		//Matrix S=dW;
		
		Matrix S= MathUtils<double>::StressVectorToTensor(rPK2_StressVector);

		double J = MathUtils<double>::Det3( rF ); // computing determinant of F

		boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
		boost::numeric::ublas::bounded_matrix<double,3,3> msaux;

		noalias(mstemp) = prod(rF,S); //computing temporary matrix F*S
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux *= J;

		if (rCauchy_StressVector.size()!=6)
		rCauchy_StressVector.resize(6);
		
		rCauchy_StressVector[0]=msaux(0,0);
		rCauchy_StressVector[1]=msaux(1,1);
		rCauchy_StressVector[2]=msaux(2,2);
		rCauchy_StressVector[3]=msaux(1,2);
		rCauchy_StressVector[4]=msaux(1,3);
		rCauchy_StressVector[5]=msaux(2,3);
	}
	
	
         int Hyperelastic3D::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo)
        {  
            if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
	    
	    if(MU.Key() == 0 || props[MU]<0.00)
                KRATOS_ERROR(std::invalid_argument,"MU has Key zero or invalid value ","");
	    
	    if(BULK_MODULUS.Key() == 0 || props[BULK_MODULUS]< 0.00) 
               KRATOS_ERROR(std::invalid_argument,"BULK_MODULUS has Key zero or invalid value ","");
    
	    return 0;
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


