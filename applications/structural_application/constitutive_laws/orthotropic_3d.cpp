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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 17:14:12 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/orthotropic_3d.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
    namespace Orthotropic3DAuxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
#ifdef _OPENMP
#pragma omp threadprivate(mstemp)
#endif
        boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
#ifdef _OPENMP
#pragma omp threadprivate(msaux)
#endif
    }
    using namespace Orthotropic3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	Orthotropic3D::Orthotropic3D() 
	: ConstitutiveLaw< Node<3> >()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Orthotropic3D::~Orthotropic3D()
	{
	}
	
	bool Orthotropic3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Orthotropic3D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
	}
	
	bool Orthotropic3D::Has( const Variable<Matrix>& rThisVariable )
	{
		if( rThisVariable == MATERIAL_DIRECTION )
			return true;
		return false;
	}
	
	double Orthotropic3D::GetValue( const Variable<double>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
	}
	
	Vector Orthotropic3D::GetValue( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
                if( rThisVariable == INTERNAL_VARIABLES )
                {
                    Vector dummy(1);
                    dummy[0] = 0.0;
                    return( dummy );
                }
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Orthotropic3D::GetValue( const Variable<Matrix>& rThisVariable )
	{
		if( rThisVariable == MATERIAL_DIRECTION )
			return mMaterialDirection;
	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
	}

	void Orthotropic3D::SetValue( const Variable<double>& rThisVariable, const double rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Orthotropic3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Orthotropic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void Orthotropic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == MATERIAL_DIRECTION )
			mMaterialDirection = rValue;
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Orthotropic3D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
        mCurrentStress = ZeroVector(6);
// 		mE = props[YOUNG_MODULUS];
// 		mNU = props[POISSON_RATIO];
		mCtangent = ZeroMatrix(6,6);
		mInSituStress = ZeroVector(6);
		mMaterialDirection = ZeroMatrix(3,3);
		mT = ZeroMatrix(6,6);

        CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);
//                 CalculateElasticMatrix(mCtangent, props[MATERIAL_PARAMETERS][0], props[MATERIAL_PARAMETERS][1]);
	}
	
	void Orthotropic3D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void Orthotropic3D::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			mInSituStress -= mCurrentStress;
			//SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
		}
	}
    
    void Orthotropic3D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {
        CalculateElasticMatrix(mCtangent, props[ORTHOTROPIC_YOUNG_MODULUS], props[ORTHOTROPIC_POISSON_RATIO], props[ORTHOTROPIC_SHEAR_RATIO]);    
    }
		
	/**
	 *	TO BE TESTED!!!
	 */
	void Orthotropic3D::CalculateElasticMatrix(Matrix& C, const array_1d<double,3>r& E, const Matrix& NU, const array_1d<double,3>r& G)
	{ 
		//setting up material matrix
		double nu_xy = NU(0,1);
		double nu_yx = NU(1,0);
		double nu_yz = NU(1,2);
		double nu_zy = NU(2,1);
		double nu_xz = NU(0,2);
		double nu_zx = NU(2,0);

		double xi = 1.00 - (nu_xy * nu_yx + nu_yz * nu_zy + nu_zx * nu_xz) - (nu_xy * nu_yz * nu_zx + nu_yx * nu_zy * nu_xz);

		double dxx = E[0] * (1.00 - nu_yz * nu_zy) / xi;
		double dxy = E[0] * (nu_xy + nu_xz * nu_zy) / xi;
		double dxz = E[0] * (nu_xz + nu_yz * nu_xy) / xi;
		double dyy = E[1] * (1.00 - nu_xz * nu_zx) / xi;
		double dyz = E[1] * (nu_yz + nu_yx * nu_xz) / xi;
		double dzz = E[2] * (1.00 - nu_yx * nu_xy) / xi;

		//filling material matrix
		C(0,0) = dxx;    C(0,1) = dxy;    C(0,2) = dxz;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
		C(1,0) = dxy;    C(1,1) = dyy;    C(1,2) = dyz;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
		C(2,0) = dxz;    C(2,1) = dyz;    C(2,2) = dzz;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
		C(3,0) = 0.0;    C(3,1) = 0.0;    C(3,2) = 0.0;    C(3,3) = G[0];  C(3,4) = 0.0;   C(3,5) = 0.0;
		C(4,0) = 0.0;    C(4,1) = 0.0;    C(4,2) = 0.0;    C(4,3) = 0.0;   C(4,4) = G[1];  C(4,5) = 0.0;
		C(5,0) = 0.0;    C(5,1) = 0.0;    C(5,2) = 0.0;    C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = G[2];

		Matrix t_matrix(6,6);

		CalculateTransformationMatrix(t_matrix);

		Matrix temp(6,6) = prod(C,t_matrix);
		noalias(C) = prod(trans(t_matrix), temp);
		
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Orthotropic3D::CalculateTransformationMatrix(Matrix& T)
	{
		Matrix& a = mMaterialDirection;

		double a11 = mMaterialDirection(0,0);
		double a12 = mMaterialDirection(0,0);
		double a13 = mMaterialDirection(0,0);

		double a21 = mMaterialDirection(1,0);
		double a22 = mMaterialDirection(1,1);
		double a23 = mMaterialDirection(1,2);

		double a31 = mMaterialDirection(2,0);
		double a32 = mMaterialDirection(2,1);
		double a33 = mMaterialDirection(2,2);

		//filling transformation matrix
		T(0,0) = a11*a11;    T(0,1) = a21*a21;    T(0,2) = a31*a31;    T(0,3) = 2.00 * a22 * a12;   	T(0,4) = 2.00 * a32 * a22;   	T(0,5) = 2.00 * a32 * a12;
		T(1,0) = a12*a12;    T(1,1) = a22*a22;    T(1,2) = a32*a32;    T(1,3) = 2.00 * a22 * a12;   	T(1,4) = 2.00 * a32 * a22;   	T(1,5) = 2.00 * a32 * a12;
		T(2,0) = a13*a13;    T(2,1) = a23*a23;    T(2,2) = a33*a33;    T(2,3) = 2.00 * a23 * a13;   	T(2,4) = 2.00 * a33 * a23;   	T(2,5) = 2.00 * a33 * a13;
		T(3,0) = a11*a12;    T(3,1) = a21*a22;    T(3,2) = a31*a32;    T(3,3) = a21 * a12 + a11 * a22;  T(3,4) = a31 * a22 + a21 * a32; T(3,5) = a31 * a12 + a11 * a32;
		T(4,0) = a12*a13;    T(4,1) = a22*a23;    T(4,2) = a32*a32;    T(4,3) = a22 * a13 + a12 * a23;  T(4,4) = a32 * a23 + a22 * a33; T(4,5) = a32 * a13 + a12 * a33;
		T(5,0) = a13*a11;    T(5,1) = a23*a21;    T(5,2) = a33*a31;    T(5,3) = a21 * a13 + a11 * a23;  T(5,4) = a31 * a23 + a21 * a33; T(5,5) = a31 * a13 + a11 * a33;
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Orthotropic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6,false);
		}
		noalias(StressVector) = prod(mCtangent,StrainVector);
		mCurrentStress = StressVector;
		noalias(StressVector) -= mInSituStress;
// 		double c1 = mE / ((1.00+mNU)*(1-2*mNU));
// 		double c2 = c1 * (1-mNU);
// 		double c3 = c1 * mNU;
// 		double c4 = c1 * 0.5 * (1 - 2*mNU);
// 
// 		StressVector[0] = c2*StrainVector[0] + c3 * (StrainVector[1] + StrainVector[2])	;
// 		StressVector[1] = c2*StrainVector[1] + c3 * (StrainVector[0] + StrainVector[2])	;
// 		StressVector[2] = c2*StrainVector[2] + c3 * (StrainVector[0] + StrainVector[1])	;
// 		StressVector[3] = c4*StrainVector[3];
// 		StressVector[4] = c4*StrainVector[4];
// 		StressVector[5] = c4*StrainVector[5];
	}
	
	/**
	 *	TO BE REVIEWED!!!
	 */
	void Orthotropic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		rResult = mCtangent;
	}
	
    //**********************************************************************
    void Orthotropic3D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

        double J = MathUtils<double>::Det3( rF );

		noalias(mstemp) = prod(rF,S);
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux *= J;

		if(rCauchy_StressVector.size() != 6)
			rCauchy_StressVector.resize(6);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(2,2);
		rCauchy_StressVector[3] = msaux(1,2);
		rCauchy_StressVector[4] = msaux(1,3);
		rCauchy_StressVector[5] = msaux(2,3);
    }
} // Namespace Kratos
