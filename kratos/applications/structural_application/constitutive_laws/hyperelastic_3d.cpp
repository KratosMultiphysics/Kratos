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
*   Date:                $Date: 2009-01-23 15:33:12 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes 
#include <iostream>
#include <limits>
// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/hyperelastic_3d.h"
#include "includes/constitutive_law.h" //nelson
#include "custom_utilities/tensor_utils.h" // nelson
#include "custom_utilities/sd_math_utils.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"


namespace Kratos
{
	// Defines a hyperelastic constitutive law in 3D space, 
	// for incompressible isotropic hyperelastic materials.
	//////////////////

    namespace Hyperelastic3DAuxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
        #pragma omp threadprivate(mstemp)
        boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
        #pragma omp threadprivate(msaux)
    }
    using namespace Hyperelastic3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	Hyperelastic3D::Hyperelastic3D() 
	: ConstitutiveLaw< Node<3> >()
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
		if( rThisVariable == INSITU_STRESS ) // I think I need to change it to PK2_STRESS_TENSOR
			return true;
		return false;
	}
	
	bool Hyperelastic3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double Hyperelastic3D::GetValue( const Variable<double>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
	}
	
	Vector Hyperelastic3D::GetValue( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Hyperelastic3D::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
	}

	void Hyperelastic3D::SetValue( const Variable<double>& rThisVariable, const double rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hyperelastic3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hyperelastic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void Hyperelastic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo ) // rCurrentProcessInfo gives info about time and Dt
 	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Hyperelastic3D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
        //mCurrentStress = ZeroVector(6);
		mMU = props[MU]; // Shear modulus
		mK = props[BULK_MODULUS]; // Bulk modulus
// 		mE = props[YOUNG_MODULUS];
// 		mNU = props[POISSON_RATIO];
		//KRATOS_WATCH(mMU);
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
    
    void Hyperelastic3D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {
		//KRATOS_TRY
        
		//CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);    
		//if (CurrentProcessInfo[TIME]/CurrentProcessInfo[DELTA_TIME] == 1.0 ) 
		//{
		CalculateConstitutiveMatrix(StrainVector, CC);  
		//}
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
	void Hyperelastic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6,false); // false para no hacer copias
		}
		//noalias(StressVector) = prod(mCtangent,StrainVector);
		//mCurrentStress = StressVector;
		//noalias(StressVector) -= mInSituStress;

			Matrix I(3,3); // creating I matrix (Identity)
			//I=new IdentityMatrix(3);
			noalias(I)=ZeroMatrix(3,3);
			I(0,0)=1.0;
			I(1,1)=1.0;
			I(2,2)=1.0;
			
				
			Matrix C(3,3); // creating matrix C (C=Ft.F)
			noalias(C)=ZeroMatrix(3,3);
			

		
				 //Version 10Feb2009
					////	////transforming "Strain Vector" into a symmetric 2nd order tensor "Voigt Notation"
					//Matrix E(3,3);
					////*Vector e(6);
					////	//noalias(e)=StrainVector;*/

					//	//SD_MathUtils<double>::VectorToTensor(StrainVector,E);
					//	E=MathUtils<double>::StressVectorToTensor(StrainVector); // Rossi's Vector form
					//	//	MathUtils<double>::VectorToTensor(Vector& StrainVector,Matrix& E);

					//	  //transforming StrainVector into Matrix form
					//	Matrix auxE(3,3);
					//	noalias(auxE)=ZeroMatrix(3,3);
					//	//Matrix auxE=MathUtils<double>::StressVectorToTensor(StrainVector);
					//	noalias(auxE)=E;
					//	auxE *= 2.0;
					//	C += auxE;
					//	C += I; //NOW I HAVE "C" ( C = 2E+I )

				//Version 10Feb2009

// put the test Strain Vector [1 0 0] and [1.0001 0 0]  /// TEST !!!!!!!!! ///
		
			//double epsilon=std::numeric_limits<double>::epsilon();
			//KRATOS_WATCH(epsilon);
			////C(0,0) = 2.0 * (1.0000+epsilon) + 1.0;
			//C(0,0) = 2.0 * 0.0000 + 1.0;
			//C(1,1) = 2.0 * 0.0000 + 1.0;
			////C(1,1) = 2.0 *(1.0000+epsilon) + 1.0;
			//C(2,2) = 2.0 * 0.0000 + 1.0;
			////C(0,1)=0.0;
			//C(0,1) = 1.0000;
			//C(1,2)=0.0;
			//C(0,2)=0.0; 
			//C(1,0)=C(0,1);
			//C(2,1)=C(1,2);
			//C(2,0)=C(0,2);

			/// END TEST !!!!!!!!! ///

			//// Version 11Feb2009
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
			//// END Version 11Feb2009

			// The following is commented as I don't have rF as input:
			//double J = MathUtils<double>::Det3( rF ); // computing determinant of F
			//noalias(mstemp) = prod(trans(rF),rF); //computing C
		
			Matrix Cinv(3,3); 
			noalias(Cinv)=ZeroMatrix(3,3);
			double det=0.0; // det of C  ( det(C)=J2 )
			//MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed
			//MathUtils<double>::InvertMatrix(C,Cinv); 
			//KRATOS_WATCH(Cinv);
			//Matrix Test(3,3); 
			//noalias(Test)=prod(C,Cinv);
			//KRATOS_WATCH(Test);

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
			W=0.5*mMU*((1.0/(pow(J2,0.3333333333333333333333333333)))*trC-3.0);
			Matrix mdwtemp(3,3);
			noalias(mdwtemp)=ZeroMatrix(3,3);
			noalias(mdwtemp)=Cinv;
			mdwtemp*=trC;
			mdwtemp*=(-1.0/3.0);
			//mdwtemp=(-1.0/3.0)*trC*Cinv;
			//Matrix I(3,3); // creating I matrix (Identity)
			//noalias(I)=ZeroMatrix(3,3);
			//I(0,0)=1.0;
			//I(1,1)=1.0;
			//I(2,2)=1.0;
			mdwtemp += I;	// parentesis solved  
			mdwtemp *= mMU;
			//mdwtemp *= (J^(2/3));
			//mdwtemp *= pow(J2,0.3333333333333333333333333333);
			mdwtemp *= 1.0/pow(J2,0.3333333333333333333333333333);
			noalias(dW) = mdwtemp; // S computed!!! :)
			// Second Piola Kirchhoff is 3x3
			
			// puting second Piola Kirchhodd in Vector form
			Vector SPKir(6); // 6 components
			SD_MathUtils<double>::TensorToVector(dW,SPKir); // Rossi's Vector form
			// Uncomment below if not considering volumetric part
			//noalias(StressVector)=SPKir;
			//noalias(rPK2_StressVector)=noalias(SPKir); // do I need this step?
			
			//***************************

			//// treating incompressibility !!!!!!!!!!!!!!!!!!!!
			
				// Introducing Ogden model pg 244 - making it quasi-incompressible
				// volumetric Strain energy

				
				Matrix Svol(3,3);
				noalias(Svol)=ZeroMatrix(3,3);
				//double k=1000000.0; //volumetric module (in Pa)
				double beta=9.0; // previous value: Beta=9
				double auxSvol=0.0;
				//double auxSvol2=0.0;
				//auxSvol=(mK/beta); // uncomment
				auxSvol=(mK/beta)*(1.0-(1.0/pow(J2,4.5))); // previous one! 
				//auxSvol-=(1.0/pow(J,beta))*(mK/beta); //CHECK THIS LINE! there is a problem here!
				//auxSvol+=(mK/beta)*(1.0-(1.0/pow(J,beta)));// uncomment
				//KRATOS_WATCH(Cinv);
				noalias(Svol)=Cinv;
				Svol*=auxSvol;
				//Svol*=2;
				//Svol=(auxSvol-auxSvol2)*Cinv;
				//KRATOS_WATCH(Svol);
				Vector vSvol(6); // to put Svol in Vector form
				SD_MathUtils<double>::TensorToVector(Svol,vSvol); // Rossi's Vector form
				// Comment if not considering volumetric part
				//KRATOS_WATCH(SPKir);
				//KRATOS_WATCH(vSvol);
				noalias(StressVector)=SPKir+vSvol; // i just changed the sign from + to -, because of -JCp (thesis, pg 21)
				//KRATOS_WATCH(StressVector);
				
				///////////////////////////////////////////NEW MATRIX///////////
				/////// TANGENT MATRIX ///////////
			Vector S1(6);
			Vector S2(6);
			if((c.size1()!=6) || c.size2() !=6)
			{
				c.resize(6,6);
			}
				if((Ctang.size1()!=6) || Ctang.size2() !=6)
			{
				Ctang.resize(6,6);
			}
			//noalias(c)=ZeroMatrix(6,6);
			CalculateStressVector (StrainVector, S1);
			//KRATOS_WATCH(S1);
			//private CalculateStressVector (const Vector& StrainVector, S1)
			int i;
				for (i=0; i<6; i++)
				{
					Vector Ep(StrainVector);
					double epsilon=std::numeric_limits<double>::epsilon();
					Ep[i]+=epsilon;
					//KRATOS_WATCH(epsilon);
					CalculateStressVector (Ep, S2);
					//KRATOS_WATCH(S2);
					int j;
					for  (j=0; j<6; j++)
					{				
						c(j,i)=(1/epsilon)*(S2[j]-S1[i]);
					}
					
				}
				noalias(Ctang)=c;
				//KRATOS_WATCH(c);
			///////////////////////////
				
			///////////// MARINO TEST //////////////
			//	// Check Derivatives //
			//	
			//double W1;
			//double W2;
			//Vector S1(6);
			//Vector S2(6);
			//Vector IntForce(6);
			//Vector IntForce_S(6);
			////Matrix c(6,6);
			////noalias(c)=ZeroMatrix(6,6);
			//CalculateEnergy (StrainVector, W1);
			//CalculateStressVector (StrainVector, S1);
			////KRATOS_WATCH(S1);
			////private CalculateStressVector (const Vector& StrainVector, S1)
			//int i;
			//	for (i=0; i<6; i++)
			//	{
			//		// pertub one degree of freedom
			//		Vector Epert(StrainVector);
			//		double pert=1e-9;
			//		Epert[i]+=pert;
			//		//KRATOS_WATCH(pert);
			//		// compute the perturbed energy and forces
			//		CalculateEnergy (Epert, W2);
			//		CalculateStressVector (Epert, S2);
			//		//KRATOS_WATCH(S2);
			//		int j;
			//		for  (j=0; j<6; j++)
			//		{				
			//			IntForce(j,0)=(1/pert)*(W2[j]-W1);
			//			//KRATOS_WATCH(IntForce);
			//			IntForce_S(j,i)=(1/pert)*(S2[j]-S1[i]); // it's not complete!
			//			// grad_Ee=grad_Ee + S*B0 * undeformed1.dvol0(ig); 
			//			//KRATOS_WATCH(IntForce_S);
			//		}
			//	}
			//	/////////////////////////////////////////
				
				//Vector InitialStress(6);
				////InitialStress[0]=36838.1;
				//InitialStress[0]=110221.0;
				////InitialStress[1]=110221.0;
				//InitialStress[1]=36838.1;
				//InitialStress[2]=110221.0;
				////InitialStress[0]=36838.1;
				//InitialStress[3]=0.0;
				//InitialStress[4]=0.0;
				//InitialStress[5]=0.0;
				////InitialStress=(36838.1, 110221.0, 110221.0, 0.0, 0.0, 0.0);
				//Vector DeltaStress(6);
				//DeltaStress=StressVector - InitialStress;
				//Vector Cline(6);
				//Cline+=DeltaStress;
				//Cline*=(1/epsilon);
				//KRATOS_WATCH(Cline);
				//double deltaE=0.0;
				//deltaE+=1+epsilon;
				//KRATOS_WATCH(deltaE);
				//double detaStress=0.0;
				//detaStress=StressVector[0];
				//KRATOS_WATCH(detaStress);
				///////////////////////////////////////////NEW MATRIX///////////

			//// computing incompressibility !!!!!!!!!!!!!!!!!!!!

			// using mixed elements:

			//	Vector EffectiveStress(6); // 6 components
			//	//EffectiveStress=J*Cinv;
			//	EffectiveStress*= J;
			//	EffectiveStress*= -1;
			//	EffectiveStress*= p; // need to define somewhere the value of p (Lagrange muliplier of the incompressibility equality constrain)
			//	//
			//***************************

			//mCurrentStress=StressVector; // What is this for?!?!?!? mCurrentStress???
			
			// Need to transform S in Sigma! How?!?! I need F !!! What the hell!!	
			//noalias(StressVector)=noalias(dW);
	}
		// PS NOALIAS SOLO EN LADO ESQUIERDO!
//*************************************************************************************************************************
		//void Hyperelastic3D::CalculateEnergy(const Vector& StrainVector, Vector& Energy)
		//{	
	
		//	Matrix I(3,3); // creating I matrix (Identity)
		//	//I=new IdentityMatrix(3);
		//	noalias(I)=ZeroMatrix(3,3);
		//	I(0,0)=1.0;
		//	I(1,1)=1.0;
		//	I(2,2)=1.0;
		//			
		//	Matrix C(3,3); // creating matrix C (C=Ft.F)
		//	noalias(C)=ZeroMatrix(3,3);
	
		//	//// IN TOTAL LAGRANGIAN! (BELOW)
		//	C(0,0) = 2.0 * StrainVector[0] + 1.0;
		//	C(1,1) = 2.0 * StrainVector[1] + 1.0;
		//	C(2,2) = 2.0 * StrainVector[2] + 1.0;
		//	C(0,1)=StrainVector[3];
		//	C(1,2)=StrainVector[4];
		//	C(0,2)=StrainVector[5]; 
		//	C(1,0)=C(0,1);
		//	C(2,1)=C(1,2);
		//	C(2,0)=C(0,2);
		//	//KRATOS_WATCH(C)
		//
		//	Matrix Cinv(3,3); 
		//	noalias(Cinv)=ZeroMatrix(3,3);
		//	double det=0.0; // det of C  ( det(C)=J2 )
		//	MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed

		//	double J2=MathUtils<double>::Det3( C ); // determinant of C is computed
		//	//KRATOS_WATCH(J2);
		//	double J=0.0;
		//	J=sqrt(J2); // J is computed
		//	//KRATOS_WATCH(J);

		//	// computing the derivative of strain energy: dW = S (the deviatoric part)
		//	
		//	double W=0.0;
		//	// Vector dW(3);
		//	Matrix dW(3,3);
		//	noalias(dW)=ZeroMatrix(3,3);
		//	double trC=0.0;
		//	trC=C(0,0)+C(1,1)+C(2,2); // computing Ckk (Trace of C = Sum of the diagonal terms) = I1 (First Invariant)
		//	W=0.5*mMU*((1.0/(pow(J2,0.3333333333333333333333333333)))*trC-3.0);
		//	
		//}

//*************************************************************************************************************************
		void Hyperelastic3D::CalculateStressVector(const Vector& StrainVector, Vector& StressVector)
		{	
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
	void Hyperelastic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)

		{
			if( rResult.size1() != 6 )
			{
			rResult.resize(6,6,false);
			}
	
	// auxiliar computations
			

			Matrix I(3,3); // creating I matrix (Identity)
			//I=new IdentityMatrix(3);
			noalias(I)=ZeroMatrix(3,3);
			I(0,0)=1.0;
			I(1,1)=1.0;
			I(2,2)=1.0;
			
			Matrix C(3,3); // creating matrix C (C=Ft.F)
			noalias(C)=ZeroMatrix(3,3);

			// Version 11Feb2009
			// IN TOTAL LAGRANGIAN! (BELOW)
			C(0,0) = 2.0 * StrainVector[0] + 1.0;
			C(1,1) = 2.0 * StrainVector[1] + 1.0;
			C(2,2) = 2.0 * StrainVector[2] + 1.0;
			C(0,1) = StrainVector[3];
			C(1,2) = StrainVector[4];
			C(0,2) = StrainVector[5]; 
			C(1,0) = C(0,1);
			C(2,1) = C(1,2);
			C(2,0) = C(0,2);

			// END Version 11Feb2009


			// The following is commented as I don't have rF as input:
			//double J = MathUtils<double>::Det3( rF ); // computing determinant of F
			//noalias(mstemp) = prod(trans(rF),rF); //computing C
		
			Matrix Cinv(3,3);
			noalias(Cinv)=ZeroMatrix(3,3);
			double det=0.0; // det of C  ( det(C)=J2 )
			//MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed
			//MathUtils<double>::InvertMatrix(C,Cinv); 

			double J2=MathUtils<double>::Det3( C ); // determinant of C is computed

	// computing the Lagrangian (or material) elasticity tensor: CC

			Matrix CC(6,6);
			noalias(CC)=ZeroMatrix(6,6); // to check matrix size!
			
			// Constant
			double J23=0.0;
			//J23=1/(J^(2/3));
			J23+=1.0/pow(J2,0.3333333333333333333333333333); 
			double constC=0.0; // 1st term. For the 2nd term: multiply by (-1/3). For the 3rd term: multiply by (-1/2)
			constC-=(1.0/3.0)*mMU*J23;
			// Volumetric part of CC
			double constV2=0.0; // Constant related to the volumetric part of CC, 2nd term invC.invC
			double constV3=0.0; // Constant related to the volumetric part of CC, 3rd term
			//double kappa;
			//kappa=1000000.0; // bulk modulus in Pa (put a high value to simulate quasi-incompressibility)
			//constV2=0.0;
			constV2-=1.0*(mK/9.0)*((11.0/2.0)*(1.0/pow(J2,1.25))-(1.0/pow(J2,4.5))); // -1 is to recover the + sign when it's summed to constC/3
			//constV2/=2.0;
			//constV2*=0.982843484;
			//KRATOS_WATCH(constV2);
			constV3+=(mK/18.0)*(1.0-(1.0/pow(J2,4.5))); // it's already considering the 1/2 related to CC+CC (3rd term), - is already included outside the parenthesis
			//constV3/=2.0;
			//constV3*=0.982843484;
			//KRATOS_WATCH(constV3);

			
			// To optimize the computation, it could be used something like the following:
					//****************************************************************************
						//CC(0,0) = constC*1constC(0)-((constC/3)+constV2)*2constC(0)+constC*;  
						
						//Matrix 1auxC=ZeroMatrix(3,3); // for (c-1*c-1), the first term
						//Matrix 2auxC=ZeroMatrix(3,3); // for (c-1*c-1+c-1*c-1), the second term
						//noalias(1constC)=prod(I,Cinv);	// is this correct?!
						//noalias(2constC)=prod(Cinv,Cinv); // is this correct?!	
					//****************************************************************************

			// this matrix is used as the elastic matrix...to compute the Cauchy Stress! (?)

				//	//filling the elasticity tensor CC (6X6) 
		

	/*	double adj=( mK /1000000.0);

		CC(0,0)= 1.00032e+006;
		CC(0,1)= 1e+006;
		CC(0,2)= 1e+006;
		CC(0,3)= 0.0;
		CC(0,4)= 0.0;
		CC(0,5)= 105.625;
		
		CC(1,0)= 1e+006;
		CC(1,1)= 1.00032e+006;
		CC(1,2)= 1e+006;
		CC(1,3)= 0.0;
		CC(1,4)= 0.0;
		CC(1,5)= 105.625;
				
		CC(2,0)= 1e+006;
		CC(2,1)= 1e+006;
		CC(2,2)= 1.00032e+006;
		CC(2,3)= 0.0;
		CC(2,4)= 0.0;
		CC(2,5)= 105.625;
				
		CC(3,0)= 0.0;
		CC(3,1)= 0.0;
		CC(3,2)= 0.0;
		CC(3,3)= 211.25;
		CC(3,4)= 2.34535e-014;
		CC(3,5)= 105.625;
				
		CC(4,0)= 0.0;
		CC(4,1)= 0.0;
		CC(4,2)= 0.0;
		CC(4,3)= 2.34535e-014;
		CC(4,4)= 211.25;
		CC(4,5)= 105.625;
				
		CC(5,0)= -105.625;
		CC(5,1)= -105.625;
		CC(5,2)= -105.625;
		CC(5,3)= -105.625;
		CC(5,4)= -105.625;
		CC(5,5)= 211.25;*/
///////////////////////////////////////
		/*CC(0,0)= Ctang(0,0);
		CC(0,1)= Ctang(0,1);
		CC(0,2)= Ctang(0,2);
		CC(0,3)= Ctang(0,3);
		CC(0,4)= Ctang(0,4);
		CC(0,5)= Ctang(0,5);
		
		CC(1,0)= Ctang(1,0);
		CC(1,1)= Ctang(1,1);
		CC(1,2)= Ctang(1,2);
		CC(1,3)= Ctang(1,3);
		CC(1,4)= Ctang(1,4);
		CC(1,5)= Ctang(1,5);
				
		CC(2,0)= Ctang(2,0);
		CC(2,1)= Ctang(2,1);
		CC(2,2)= Ctang(2,2);
		CC(2,3)= Ctang(2,3);
		CC(2,4)= Ctang(2,4);
		CC(2,5)= Ctang(2,5);
				
		CC(3,0)= Ctang(3,0);
		CC(3,1)= Ctang(3,1);
		CC(3,2)= Ctang(3,2);
		CC(3,3)= Ctang(3,3);
		CC(3,4)= Ctang(3,4);
		CC(3,5)= Ctang(3,5);
				
		CC(4,0)= Ctang(4,0);
		CC(4,1)= Ctang(4,1);
		CC(4,2)= Ctang(4,2);
		CC(4,3)= Ctang(4,3);
		CC(4,4)= Ctang(4,4);
		CC(4,5)= Ctang(4,5);
				
		CC(5,0)= Ctang(5,0);
		CC(5,1)= Ctang(5,1);
		CC(5,2)= Ctang(5,2);
		CC(5,3)= Ctang(5,3);
		CC(5,4)= Ctang(5,4);
		CC(5,5)= Ctang(5,5);*/


		//noalias(rResult)=CC*adj;
		
		//noalias(rResult)=CC;
		/////////// TANGENT MATRIX ///////////
		//	Vector S1(6);
		//	Vector S2(6);
		//	if((c.size1()!=6) || c.size2() !=6)
		//	{
		//		c.resize(6,6);
		//	}
		//		if((Ctang.size1()!=6) || Ctang.size2() !=6)
		//	{
		//		Ctang.resize(6,6);
		//	}
		//	//noalias(c)=ZeroMatrix(6,6);
		//	CalculateStressVector (StrainVector, S1);
		//	//KRATOS_WATCH(S1);
		//	//private CalculateStressVector (const Vector& StrainVector, S1)
		//	int i;
		//		for (i=0; i<6; i++)
		//		{
		//			Vector Ep(StrainVector);
		//			double epsilon=std::numeric_limits<double>::epsilon();
		//			Ep[i]+=epsilon;
		//			//KRATOS_WATCH(epsilon);
		//			CalculateStressVector (Ep, S2);
		//			//KRATOS_WATCH(S2);
		//			int j;
		//			for  (j=0; j<6; j++)
		//			{				
		//				c(j,i)=(1/epsilon)*(S2[j]-S1[i]);
		//			}
		//			
		//		}
		//		noalias(Ctang)=c;
		//noalias(rResult)=Ctang;
		////KRATOS_WATCH (rResult);


//	 Rossi notation, with right constants

			// 1ST ROW  ----------> refference xx
			//CC(0,0) = constC*1constC(0)-((constC/3)+constV2)*2constC(0)+constC*;  
			 CC(0,0) = constC*I(0,0)*Cinv(0,0)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(0,0)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(0,0)+Cinv(0,0)*Cinv(0,0));
				CC(0,1) = constC*I(0,0)*Cinv(1,1)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(1,1)-((constC/2.0)+constV3)*(Cinv(0,1)*Cinv(0,1)+Cinv(0,1)*Cinv(0,1));  
					CC(0,2) = constC*I(0,0)*Cinv(2,2)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(2,2)-((constC/2.0)+constV3)*(Cinv(0,2)*Cinv(0,2)+Cinv(0,2)*Cinv(0,2));    
						CC(0,3) = constC*I(0,0)*Cinv(0,1)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(0,1)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(0,1)+Cinv(0,1)*Cinv(0,0));   
							CC(0,4) = constC*I(0,0)*Cinv(1,2)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(1,2)-((constC/2.0)+constV3)*(Cinv(0,1)*Cinv(0,2)+Cinv(0,2)*Cinv(0,1));   
								CC(0,5) = constC*I(0,0)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(0,0)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(0,2)+Cinv(0,2)*Cinv(0,0));   
			// 2ND ROW  ----------> refference yy
				CC(1,0) = CC(0,1);   
					CC(1,1) = constC*I(1,1)*Cinv(1,1)-((constC/3.0)+constV2)*Cinv(1,1)*Cinv(1,1)-((constC/2.0)+constV3)*(Cinv(1,1)*Cinv(1,1)+Cinv(1,1)*Cinv(1,1));
						 CC(1,2) = constC*I(1,1)*Cinv(2,2)-((constC/3.0)+constV2)*Cinv(1,1)*Cinv(2,2)-((constC/2.0)+constV3)*(Cinv(1,2)*Cinv(1,2)+Cinv(1,2)*Cinv(1,2));   
							CC(1,3) = constC*I(1,1)*Cinv(0,1)-((constC/3.0)+constV2)*Cinv(1,1)*Cinv(0,1)-((constC/2.0)+constV3)*(Cinv(1,0)*Cinv(1,1)+Cinv(1,1)*Cinv(1,0));   
								CC(1,4) = constC*I(1,1)*Cinv(1,2)-((constC/3.0)+constV2)*Cinv(1,1)*Cinv(1,2)-((constC/2.0)+constV3)*(Cinv(1,1)*Cinv(1,2)+Cinv(1,2)*Cinv(1,1));   
									CC(1,5) = constC*I(1,1)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(1,1)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(1,0)*Cinv(1,2)+Cinv(1,2)*Cinv(1,0));   
			// 3RD ROW  ----------> refference zz
			CC(2,0) = CC(0,2);		
				CC(2,1) = CC(1,2);    
					CC(2,2) = constC*I(2,2)*Cinv(2,2)-((constC/3.0)+constV2)*Cinv(2,2)*Cinv(2,2)-((constC/2.0)+constV3)*(Cinv(2,2)*Cinv(2,2)+Cinv(2,2)*Cinv(2,2));    
						CC(2,3) = constC*I(2,2)*Cinv(0,1)-((constC/3.0)+constV2)*Cinv(2,2)*Cinv(0,1)-((constC/2.0)+constV3)*(Cinv(2,0)*Cinv(2,1)+Cinv(2,1)*Cinv(2,0)); 
							CC(2,4) = constC*I(2,2)*Cinv(1,2)-((constC/3.0)+constV2)*Cinv(2,2)*Cinv(1,2)-((constC/2.0)+constV3)*(Cinv(2,1)*Cinv(2,2)+Cinv(2,1)*Cinv(2,2)); 
								CC(2,5) = constC*I(2,2)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(2,2)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(2,0)*Cinv(2,2)+Cinv(2,0)*Cinv(2,2)); 
			// 4TH ROW  ----------> refference to xy
			CC(3,0) = CC(0,3);   
				CC(3,1) = CC(1,3);   
					CC(3,2) = CC(2,3);      
						CC(3,3) = constC*I(0,1)*Cinv(0,1)-((constC/3.0)+constV2)*Cinv(0,1)*Cinv(0,1)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(1,1)+Cinv(0,1)*Cinv(1,0));    
							CC(3,4) = constC*I(0,1)*Cinv(1,2)-((constC/3.0)+constV2)*Cinv(0,1)*Cinv(1,2)-((constC/2.0)+constV3)*(Cinv(0,1)*Cinv(1,2)+Cinv(0,2)*Cinv(1,1));
								CC(3,5) = constC*I(0,1)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(0,1)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(1,2)+Cinv(0,2)*Cinv(1,0));
			// 5TH ROW  ----------> refference to xz
			CC(4,0) = CC(0,4);   
				CC(4,1) = CC(1,4);   
					CC(4,2) = CC(2,4);   
						CC(4,3) = CC(3,4);   
							CC(4,4) = constC*I(1,2)*Cinv(1,2)-((constC/3.0)+constV2)*Cinv(1,2)*Cinv(1,2)-((constC/2.0)+constV3)*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(2,1)); 
								CC(4,5) = constC*I(1,2)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(1,2)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(1,0)*Cinv(2,2)+Cinv(1,2)*Cinv(2,0)); 
			// 6TH ROW  ----------> refference to yz
			CC(5,0) = CC(0,5);   
				CC(5,1) = CC(1,5);   
					CC(5,2) = CC(2,5);   
						CC(5,3) = CC(3,5);   
							CC(5,4) = CC(4,5);    
								CC(5,5) = constC*I(0,2)*Cinv(0,2)-((constC/3.0)+constV2)*Cinv(0,2)*Cinv(0,2)-((constC/2.0)+constV3)*(Cinv(0,0)*Cinv(2,2)+Cinv(0,2)*Cinv(2,0)); 
 

		noalias(rResult)=CC;
		//KRATOS_WATCH(CC);

	}

	
	//*************************************************************************************************************************
	void Hyperelastic3D::CalculateCauchyStresses(	
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector, // in my case I already have S in a matrix form! How to read it?
		const Vector& rGreenLagrangeStrainVector)
	{
		//Matrix S=MathUtils<double>::StressVectorToTensor(rPK2_StressVector); // I think I don't need this!!! S is already a matrix!
		//Matrix S=dW;
		
		Matrix S= MathUtils<double>::StressVectorToTensor(rPK2_StressVector);
		//KRATOS_WATCH(S);
		double J = MathUtils<double>::Det3( rF ); // computing determinant of F
		//double J = MathUtils<double>::Det2( rF );
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
		//KRATOS_WATCH(rCauchy_StressVector);
	}



} // Namespace Kratos
