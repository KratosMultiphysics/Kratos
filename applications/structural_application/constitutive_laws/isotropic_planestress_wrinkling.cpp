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
 
//    
//   Project Name:        Kratos        
//   Last modified by:    $Author: rrossi $ 
//   Date:                $Date: 2007-11-22 11:49:23 $ 
//   Revision:            $Revision: 1.1 $ 
// 
// 
 
  
// System includes  
 
 
// External includes  
 
 
// Project includes  
#include<cmath> 
#include "includes/define.h" 
#include "constitutive_laws/isotropic_planestress_wrinkling.h" 
#include "utilities/math_utils.h" 
#include "includes/variables.h" 
#include "includes/process_info.h" 
#include "includes/properties.h" 
#include "structural_application.h" 
 
namespace Kratos 
{ 
 
	IsotropicPlaneStressWrinkling::IsotropicPlaneStressWrinkling() 
	: ConstitutiveLaw()
	{ 
		mEw.resize(3); 
	} 
 
	IsotropicPlaneStressWrinkling::IsotropicPlaneStressWrinkling(const double E, const double NU) 
	: ConstitutiveLaw()
	{  
		mEw.resize(3); 
		mE = E; 
		mNU = NU; 
 
		mCtangent = CalculateElasticMatrix(mE,mNU); 
	} 
 
	IsotropicPlaneStressWrinkling::~IsotropicPlaneStressWrinkling() 
	{ 
	} 
 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::InitializeMaterial( 	const Properties& props,
				const GeometryType& geom,
     				const Vector& ShapeFunctionsValues )
	{  
		
		mE = props[YOUNG_MODULUS]; 
		mNU = props[POISSON_RATIO]; 
 
		mCtangent = CalculateElasticMatrix(mE,mNU); 
	} 
 
	//********************************************************************** 
	Matrix IsotropicPlaneStressWrinkling::CalculateElasticMatrix(const double E, const double NU) 
	{  
		Matrix C(3,3); 
		//C.Resize(3,3); 
		double c1 = E / (1.00 - NU*NU); 
		double c2 = c1 * NU; 
		double c3 = 0.5*E / (1 + NU); 
 
		C(0,0) = c1;	C(0,1) = c2;	C(0,2) = 0.0; 
		C(1,0) = c2;	C(1,1) = c1;	C(1,2) = 0.0; 
		C(2,0) = 0.0;	C(2,1) = 0.0;	C(2,2) = c3; 
 
		return C; 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateStress( const Vector& StrainVector, 
			Vector& StressVector)
	{ 
		if(StressVector.size()!=3)
			StressVector.resize(3,false);
		noalias(StressVector) = prod(mCtangent,StrainVector); 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateConstitutiveMatrix( const Vector& StrainVector, 
			Matrix& ElasticityTensor )
	{ 
		if(ElasticityTensor.size1() != 3)
			ElasticityTensor.resize(3,3,false);
		noalias(ElasticityTensor) = mCtangent;
	} 
 
 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateStressEig(const Vector& StressVector,double& smin, Vector& mineigenvect,double& smax, Vector& maxeigenvect) 
	{ 
		const double& s11 = StressVector[0]; 
		const double& s22 = StressVector[1]; 
		const double& s12 = StressVector[2]; 
		mineigenvect.resize(2); 
		maxeigenvect.resize(2); 
		if (s12 != 0) 
		{ 
			smin = 0.5*(s22  + s11  - sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12))); 
			mineigenvect[0] = -(s22 / 2 - s11 / 2 + sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12)) / 2); 
			mineigenvect[1] = s12; 
			double norm = sqrt(mineigenvect[0]*mineigenvect[0] + mineigenvect[1]*mineigenvect[1]); 
			mineigenvect = mineigenvect/norm; 
 
			smax = 0.5*(s22  + s11  + sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12))); 
		} 
		else if (s11 > s22) 
		{ 
			smin = s22; 
			smax = s11; 
			mineigenvect[0] = 0.0; 
			mineigenvect[1] = 1.0; 
		} 
		else 
		{ 
			smin = s11; 
			smax = s22; 
			mineigenvect[0] = 1.0; 
			mineigenvect[1] = 0.0; 
		} 
 
		maxeigenvect[0] = - mineigenvect[1]; 
		maxeigenvect[1] =   mineigenvect[0]; 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateStressEigNonNormalized(const Vector& StressVector,double& smin, Vector& mineigenvect,double& smax, Vector& maxeigenvect) 
	{ 
		const double& s11 = StressVector[0]; 
		const double& s22 = StressVector[1]; 
		const double& s12 = StressVector[2]; 
		mineigenvect.resize(2); 
		maxeigenvect.resize(2); 
		if (s12 != 0) 
		{ 
			smin = 0.5*(s22  + s11  - sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12))); 
			mineigenvect[0] = -(s22 / 2 - s11 / 2 + sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12)) / 2); 
			mineigenvect[1] = s12; 
			smax = 0.5*(s22  + s11  + sqrt(s22*s22 - 2 * s11 * s22 + s11*s11 + (4 * s12*s12))); 
		} 
		else if (s11 > s22) 
		{ 
			smin = s22; 
			smax = s11; 
			mineigenvect[0] = 0.0; 
			mineigenvect[1] = 1.0; 
		} 
		else 
		{ 
			smin = s11; 
			smax = s22; 
			mineigenvect[0] = 1.0; 
			mineigenvect[1] = 0.0; 
		} 
 
		maxeigenvect[0] = - mineigenvect[1]; 
		maxeigenvect[1] =   mineigenvect[0]; 
	} 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::PrincipSTRAIN(const Vector& StrainVector,double& eps1, double& eps2) 
	{ 
		const double& e1 =  StrainVector[0]; 
		const double& e2 = StrainVector[1]; 
		const double& e3 = StrainVector[2]; 
 
		eps1 = e2 / 2 + e1 / 2 + sqrt(e2*e2 - 2 * e2 * e1 + e1*e1 + e3*e3) / 2; 
		eps2 = e2 / 2 + e1 / 2 - sqrt(e2*e2 - 2 * e2 * e1 + e1*e1 + e3*e3) / 2; 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::PrincipSTRESS(const Vector& StressVector,double& str1, double& str2) 
	{ 
		const double& s1 =  StressVector[0]; 
		const double& s2 = StressVector[1]; 
		const double& s3 = StressVector[2]; 
 
		str1 = s2 / 2 + s1 / 2 + sqrt(s2*s2 - 2 * s2 * s1 + s1*s1 + (4 * s3*s3)) / 2; 
		str2 = s2 / 2 + s1 / 2 - sqrt(s2*s2 - 2 * s2 * s1 + s1*s1 + (4 * s3*s3)) / 2;	 
	} 
 
	//********************************************************************** 
	unsigned int IsotropicPlaneStressWrinkling::AssessState(const Vector& StressVector,const Vector& StrainVector) 
	{ 
		double e1, e2; 
		double s1, s2; 
		unsigned int state; 
		PrincipSTRAIN(StrainVector,e1,e2); 
		PrincipSTRESS(StressVector,s1,s2); 
 
		if (s2>=0) 
			state = 0; //taut 
		else if (e1 >0) 
			state = 1; //wrinkled 
		else  
			state = 2; //slack 
 
		return state; 
	} 
    
    void IsotropicPlaneStressWrinkling::CalculateMaterialResponse( const Vector& StrainVector,
                                      const Matrix& DeformationGradient,
                                      Vector& StressVector,
                                      Matrix& AlgorithmicTangent,
                                      const ProcessInfo& CurrentProcessInfo,
                                      const Properties& props, 
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      bool CalculateStresses,
                                      int CalculateTangent,
                                      bool SaveInternalVariables )
    {
        CalculateStress(StrainVector, StressVector);
        CalculateConstitutiveMatrix(StrainVector, AlgorithmicTangent);
    }
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateTangentMatrix(const Vector& StrainVector) 
	{ 
		Matrix Celastic = CalculateElasticMatrix(mE,mNU);  
		Vector ElasticStress(3); 
		noalias(ElasticStress) = prod(Celastic,StrainVector); 
//		Matrix::MatVec(ElasticStress,Celastic,StrainVector); 
 
		//the state is determined once per time step if the flag = 1 otherwise once per iteration 
		unsigned int state = AssessState(ElasticStress,StrainVector); 
 
		if (state == 0) // taut 
		{ 
			mCtangent = Celastic; 
			mEw = ZeroVector(3); 
		} 
		else if (state == 1) //wrinkled 
		{ 
			double smin=0.00; double smax=0.00; 
			Vector vmin; Vector vmax; 
			CalculateStressEig(ElasticStress,smin,vmin,smax,vmax); 
			mCtangent = ConstructUnidirectionalConstitutiveMatrix(vmax); 
			mCtangent *= mE; 
 
double E11 = StrainVector[0]; 
double E22 = StrainVector[1]; 
double E12 = StrainVector[2]; 
double theta = 0.5*atan2(E12,(E11-E22)); 
Vector v(2); 
double c = cos(theta); double s = sin(theta); 
v[0] = c; v[1] = s; 
 
Matrix B(3,2); 
B(0,0) = 2*v[0]; B(0,1) = 0.00; 
B(1,0) = 0.00; B(1,1) = 2*v[1]; 
B(2,0) = v[1]; B(2,1) = v[0]; 
 
double denom = 2.00*(E11*E11-2*E11*E22+E22*E22+E12*E12); 
Matrix A(2,3); 
A(0,0) = s*E12; A(0,1) = -s*E12; A(0,2) = -s*(E11-E22); 
A(1,0) = -c*E12; A(1,1) = c*E12; A(1,2) = c*(E11-E22); 
A /= denom; 
 
Matrix temp(3,3); 
noalias(temp) = prod(B,A); 
//Matrix::MatMul(temp,B,A); 
double emax,emin; 
PrincipSTRAIN(StrainVector,emax,emin); 
temp *= mE*emax; 
mCtangent += temp; 
//std::cout << "*********************"  << std::endl; 
//std::cout << mCtangent  << std::endl; 
//std::cout << temp  << std::endl; 
 
			//Matrix ShearTerm(3,3); 
			//Vector ShearTermBase(3); 
			//ShearTermBase[0] = vmin[0]*vmax[0]; 
			//ShearTermBase[1] = vmin[1]*vmax[1]; 
			//ShearTermBase[2] = 0.5*(vmin[0]*vmax[1]+vmin[1]*vmax[0]); 
			//ShearTerm = MathUtils<double>::TensorProduct3(ShearTermBase,ShearTermBase); 
			//ShearTerm *= mE/((1+mNU)); 
 
//std::cout << ShearTerm  << std::endl; 
//std::cout << "*********************"  << std::endl; 
 
			//mCtangent += ShearTerm; 
 
			//calculation of the wrinkling strain 
			Matrix Cinv; 
			double det; 
			MathUtils<double>::InvertMatrix3(Celastic,Cinv,det); 
			Vector EffectiveStress(3); 
		noalias(EffectiveStress) = prod(mCtangent,StrainVector); 
		//	Matrix::MatVec(EffectiveStress,mCtangent,StrainVector); 
		noalias(mEw) = prod(Cinv,EffectiveStress); 
		//	Matrix::MatVec(mEw,Cinv,EffectiveStress); 
			mEw -= StrainVector; 
 
		} 
		else if (state == 2) // slack 
		{ 
			mCtangent.resize(3,3); 
			mCtangent = ZeroMatrix(3); 
 
			////shear terms 
			//double smin=0.00; double smax=0.00; 
			//Vector vmin; Vector vmax; 
			//CalculateStressEig(ElasticStress,smin,vmin,smax,vmax); 
			//Matrix ShearTerm(3,3); 
			//Vector ShearTermBase(3); 
			//ShearTermBase[0] = vmin[0]*vmax[0]; 
			//ShearTermBase[1] = vmin[1]*vmax[1]; 
			//ShearTermBase[2] = 0.5*(vmin[0]*vmax[1]+vmin[1]*vmax[0]); 
			//ShearTerm = MathUtils<double>::TensorProduct3(ShearTermBase,ShearTermBase); 
			//ShearTerm *= mE/((1+mNU)); 
			//mCtangent += ShearTerm; 
 
			mEw[0] = -StrainVector[0]; 
			mEw[1] = -StrainVector[1]; 
			mEw[2] = -StrainVector[2]; 
		} 
 
	} 
 
 
	//********************************************************************** 
	Matrix IsotropicPlaneStressWrinkling::ConstructUnidirectionalConstitutiveMatrix(const Vector& v) 
	{ 
		Matrix A(3,3); 
 
		A(0,0) = v[0]*v[0]*v[0]*v[0];      A(0,1) = v[0]*v[0]*v[1]*v[1];      A(0,2) = v[0]*v[0]*v[0]*v[1]; 
		A(1,0) = v[1]*v[1]*v[0]*v[0];      A(1,1) = v[1]*v[1]*v[1]*v[1];      A(1,2) = v[1]*v[1]*v[0]*v[1]; 
		A(2,0) = v[0]*v[1]*v[0]*v[0];      A(2,1) = v[0]*v[1]*v[1]*v[1];      A(2,2) = v[1]*v[0]*v[0]*v[1]; 
 
		return A; 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::UpdateMaterial( const Vector& StrainVector,
						const Properties& props,
						const GeometryType& geom,
						const Vector& ShapeFunctionsValues ,
						const ProcessInfo& CurrentProcessInfo)

	{ 
		KRATOS_TRY 
		
// 				KRATOS_WATCH("updating the material");
/*			int wrinkling_approach = CurrentProcessInfo[WRINKLING_APPROACH]; 
		int it = CurrentProcessInfo[NL_ITERATION_NUMBER];  */
//		int solstep = CurrentProcessInfo.GetCurrentSolutionStep(); 
 
		CalculateTangentMatrix(StrainVector); 
		
// 		if(wrinkling_approach == 1 ) 
// 		{ 
// 			if(it == 1) 
// 			{ 
// 				CalculateTangentMatrix(StrainVector); 
// 			} 
// 		} 
// 		else if(wrinkling_approach == 2 ) 
// 		{ 
// 			CalculateTangentMatrix(StrainVector); 
// 		} 
// 		else 
// 		{ 
// 			mCtangent = CalculateElasticMatrix(mE,mNU);  
// 		} 
 
//KRATOS_WATCH(mCtangent); 
		KRATOS_CATCH("") 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::Calculate( const Variable<Matrix >& rVariable, 
			Matrix& Output, 
   			const ProcessInfo& rCurrentProcessInfo)
	{ 
		if(rVariable==AUXILIARY_MATRIX_1) 
		{ 
			if(Output.size1() != 1 || Output.size2() != 3)
				Output.resize(1,3,false);
			Output(0,0) = mEw[0]; 
			Output(0,1) = mEw[1]; 
			Output(0,2) = mEw[2]*0.5; 
		} 
		KRATOS_ERROR(std::logic_error,"Trying to Calculate an inexisting variable" , "") 
	} 
 
	//********************************************************************** 
	void IsotropicPlaneStressWrinkling::CalculateCauchyStresses( 
		Vector& Cauchy_StressVector, 
		const Matrix& F, 
		const Vector& PK2_StressVector, 
		const Vector& GreenLagrangeStrainVector) 
	{  
		//compute effective thickness variation 
		Vector EffectiveStrain =  GreenLagrangeStrainVector; 
 
		Matrix Feff(2,2); 
		if (MathUtils<double>::Norm3(mEw) != 0) 
		{ 
			EffectiveStrain += mEw; 
			//double emin=0.00; double emax=0.00; 
			Vector vmin; Vector vmax; 
			Matrix Finv(2,2); 
			double Fdet; 
			MathUtils<double>::InvertMatrix2(F,Finv,Fdet); 
			Matrix wwt(2,2); 
			Matrix Ewrink(2,2); 
			Ewrink(0,0) = mEw[0];      Ewrink(0,1) = mEw[2]*0.5; 
			Ewrink(1,0) = mEw[2]*0.5;  Ewrink(1,1) = mEw[1]; 
 
			noalias(wwt) = prod( Ewrink, trans(Finv) ); 
			wwt = prod(Finv, wwt); std::cout << "oribbbbbbbbbile" << std::endl; 
//			noalias(wwt) = prod( Finv, prod( Ewrink, trans(Finv) ) ); std::cout << "efficiency" << std::endl; 
//			noalias(wwt) = prod( Finv, prod<Matrix>(Ewrink,trans(Finv)) ); 
//			Matrix::MatMulAndAdd_B_D_Btrans(wwt,Finv,Ewrink); 
			double alpha = fabs(wwt(0,0) + wwt(1,1)); 
			wwt /= alpha; 
			double beta = -1.00 + sqrt(2*alpha+1.00); 
			noalias(Feff) = prod(wwt,F); 
			//Matrix::MatMul(Feff,wwt,F); 
			Feff *= beta; 
			Feff += F; 
		} 
		else 
		{ 
			Feff = F; 
		} 
 
		double epsz = ((mNU)/(1-mNU))*(EffectiveStrain[0] + EffectiveStrain[1]); 
		double h_h0_ratio = sqrt(1-2.00*epsz); 
		double DetF = MathUtils<double>::Det2(Feff); 
		Matrix S(2,2);  
		S(0,0) = PK2_StressVector[0]; S(0,1) = PK2_StressVector[2]; 
		S(1,0) = PK2_StressVector[2]; S(1,1) = PK2_StressVector[1]; 
		Matrix s(2,2); 
 
		noalias(s) = prod(Feff, Matrix(prod(S, trans(Feff)) ) ); std::cout << "efficienza!!" << std::endl; 
//		noalias(s) = prod(Feff, prod<Matrix>(S, trans(Feff) ) ); 
//		Matrix::MatMulAndAdd_B_D_Btrans(s,Feff,S); //s = F*S*Ftrans 
		Cauchy_StressVector[0] = s(0,0); Cauchy_StressVector[1] = s(1,1); Cauchy_StressVector[2] = s(0,1); 
		Cauchy_StressVector /= DetF * h_h0_ratio; 
 
	} 
 
	//********************************************************************** 
	//returns h/h0 
	double IsotropicPlaneStressWrinkling::CalculateThicknessRatio( 
		const Vector& GreenLagrangeStrainVector) 
	{ 
		double epsz = ((mNU)/(1-mNU))*(GreenLagrangeStrainVector[0] + GreenLagrangeStrainVector[1]); 
		return sqrt(1-2.00*epsz); 
	} 
	
	int IsotropicPlaneStressWrinkling::Check(const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY
          
            
            if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
            
            if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

	    const double& nu = props[POISSON_RATIO];
	    const bool check = bool( (nu >0.999 && nu<1.01 ) || (nu < -0.999 && nu > -1.01 ) );
	    if(POISSON_RATIO.Key() == 0 || check==true) 
                KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

	    return 0;
	    
            KRATOS_CATCH("");
        }
 
 
} // Namespace Kratos 
 
 
