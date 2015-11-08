/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
*   Last Modified by:    $Author:   Massimo Petracca$
*   Date:                $Date:     19-09-2014$
*   Revision:            $Revision: 1.0$
*
* ***********************************************************/

#include "damage_tc_3d_law.h"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"
#include "includes/variables.h"

#include "custom_utilities/imperfection_utilities.h"

#include <math.h>

#define DAM_TC_PREC 1.0E-12

#define DAM_TC_SIGN(X) (X == 0.0 ? 0.0 : ( X > 0.0 ? 1.0 : -1.0 ))
#define DAM_TC_HVS(X) ( X > 0.0 ? 1.0 : ( X < 0.0 ? -1.0 : 0.5) )
#define DAM_TC_POS(X) ( X > 0.0 ? X : 0.0 )

//#define DAM_TC_TEST_RANKINE
//#define DAM_TC_SECANT

//#define DAM_TC_NO_REGULARIZATION
#ifdef DAM_TC_NO_REGULARIZATION
#define DAM_TC_FIXED_LCH 200.0
#endif // DAM_TC_NO_REGULARIZATION

//#define DAM_TC_TEST_TENSILE_LAW

#define DAM_TC_OPTIMIZE_LCH

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // !M_PI

namespace Kratos
{

	// =====================================================================================================
	//
	// UTILITIES
	//
	// =====================================================================================================

#ifdef DAM_TC_OPTIMIZE_LCH

	inline double get_optimized_lch(const Element::GeometryType& geom)
	{
		double lch = geom.Length();
		if(geom.WorkingSpaceDimension() == 3) {
			if(geom.PointsNumber() == 4 /*&& geom.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral*/) {
				// 3D 4N (shell)
				double ax = (geom[0].X0()+geom[3].X0())/2.0;
				double ay = (geom[0].Y0()+geom[3].Y0())/2.0;
				double az = (geom[0].Z0()+geom[3].Z0())/2.0;
				double bx = (geom[1].X0()+geom[2].X0())/2.0;
				double by = (geom[1].Y0()+geom[2].Y0())/2.0;
				double bz = (geom[1].Z0()+geom[2].Z0())/2.0;
				double cx = (geom[0].X0()+geom[1].X0())/2.0;
				double cy = (geom[0].Y0()+geom[1].Y0())/2.0;
				double cz = (geom[0].Z0()+geom[1].Z0())/2.0;
				double dx = (geom[2].X0()+geom[3].X0())/2.0;
				double dy = (geom[2].Y0()+geom[3].Y0())/2.0;
				double dz = (geom[2].Z0()+geom[3].Z0())/2.0;
				double v1x = bx-ax;
				double v1y = by-ay;
				double v1z = bz-az;
				double v2x = dx-cx;
				double v2y = dy-cy;
				double v2z = dz-cz;
				double lx = std::sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
				double ly = std::sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
				lch = std::min(lx,ly);
			}
		}
		return lch;
	}

#endif // DAM_TC_OPTIMIZE_LCH

	inline double b3_eval_bezier(double xi,
								 double x0, double x1, double x2,
								 double y0, double y1, double y2)
	{
		double A = x0-2.0*x1+x2;
		double B = 2.0*(x1-x0);
		double C = x0-xi;
		if(std::abs(A)<1.0E-12)
		{
			x1 = x1+1.0E-6*(x2-x0);
			A = x0-2.0*x1+x2;
			B = 2.0*(x1-x0);
			C = x0-xi;
		}
		double D = B*B-4.0*A*C;
		double t = (-B+std::sqrt(D))/(2.0*A);
		return (y0-2.0*y1+y2)*t*t+2.0*(y1-y0)*t+y0;
	}

	inline double b3_eval_area(double x1,double x2,double x3,double y1,double y2,double y3)
	{
		return x2*y1/3.0 + x3*y1/6.0 - x2*y3/3.0 + x3*y2/3.0 + x3*y3/2.0
			  -x1*(y1/2.0 + y2/3.0 + y3/6.0);
	}

	inline double b3_calc_G(double sp, double sk, double sr,
							double ep, double ej, double ek, double er, double eu)
	{
		double G1 = ep*sp/2.0;
		double G2 = b3_eval_area(ep,ej,ek,sp,sp,sk);
		double G3 = b3_eval_area(ek,er,eu,sk,sr,sr);
		return G1+G2+G3;
	}

	inline void b3_stretch(double S, double ep, double &ej, double &ek, double &er, double &eu)
	{
		ej = ej+(ej-ep)*S;
		ek = ek+(ek-ep)*S;
		er = er+(er-ep)*S;
		eu = eu+(eu-ep)*S;
	}

	inline void dam_tc_swap(double& a, double& b) {
		double c = a; a = b; b = c;
	}
	inline void dam_tc_swap(array_1d<double,3>& a, array_1d<double,3>& b, array_1d<double,3>& c) {
		noalias(c)=a; noalias(a)=b; noalias(b)=c;
	}
	inline void dam_tc_add_dyad(const array_1d<double,3>& n, double x, Vector& S) {
		S(0)+=x*n(0)*n(0);
		S(1)+=x*n(1)*n(1);
		S(2)+=x*n(2)*n(2);
		S(3)+=x*n(0)*n(1);
		S(4)+=x*n(1)*n(2);
		S(5)+=x*n(0)*n(2);
	}
	
	/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
	domain Java Matrix library JAMA. */

	#ifdef EIG3_MAX
	#undef EIG3_MAX
	#endif
	#define EIG3_MAX(a, b) ((a)>(b)?(a):(b))

	inline double hypot2(double x, double y) {
		return std::sqrt(x*x+y*y);
	}

	// Symmetric Householder reduction to tridiagonal form.
	inline void tred2(Matrix& V, Vector& d, Vector& e) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.
	  int n = d.size();
	  for (int j = 0; j < n; j++) {
		d(j) = V(n-1,j);
	  }

	  // Householder reduction to tridiagonal form.

	  for (int i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
		  scale = scale + fabs(d(k));
		}
		if (scale == 0.0) {
		  e(i) = d(i-1);
		  for (int j = 0; j < i; j++) {
			d(j) = V(i-1,j);
			V(i,j) = 0.0;
			V(j,i) = 0.0;
		  }
		} else {

		  // Generate Householder vector.

		  for (int k = 0; k < i; k++) {
			d(k) /= scale;
			h += d(k) * d(k);
		  }
		  double f = d(i-1);
		  double g = sqrt(h);
		  if (f > 0) {
			g = -g;
		  }
		  e(i) = scale * g;
		  h = h - f * g;
		  d(i-1) = f - g;
		  for (int j = 0; j < i; j++) {
			e(j) = 0.0;
		  }

		  // Apply similarity transformation to remaining columns.

		  for (int j = 0; j < i; j++) {
			f = d(j);
			V(j,i) = f;
			g = e(j) + V(j,j) * f;
			for (int k = j+1; k <= i-1; k++) {
			  g += V(k,j) * d(k);
			  e(k) += V(k,j) * f;
			}
			e(j) = g;
		  }
		  f = 0.0;
		  for (int j = 0; j < i; j++) {
			e(j) /= h;
			f += e(j) * d(j);
		  }
		  double hh = f / (h + h);
		  for (int j = 0; j < i; j++) {
			e(j) -= hh * d(j);
		  }
		  for (int j = 0; j < i; j++) {
			f = d(j);
			g = e(j);
			for (int k = j; k <= i-1; k++) {
			  V(k,j) -= (f * e(k) + g * d(k));
			}
			d(j) = V(i-1,j);
			V(i,j) = 0.0;
		  }
		}
		d(i) = h;
	  }

	  // Accumulate transformations.

	  for (int i = 0; i < n-1; i++) {
		V(n-1,i) = V(i,i);
		V(i,i) = 1.0;
		double h = d(i+1);
		if (h != 0.0) {
		  for (int k = 0; k <= i; k++) {
			d(k) = V(k,i+1) / h;
		  }
		  for (int j = 0; j <= i; j++) {
			double g = 0.0;
			for (int k = 0; k <= i; k++) {
			  g += V(k,i+1) * V(k,j);
			}
			for (int k = 0; k <= i; k++) {
			  V(k,j) -= g * d(k);
			}
		  }
		}
		for (int k = 0; k <= i; k++) {
		  V(k,i+1) = 0.0;
		}
	  }
	  for (int j = 0; j < n; j++) {
		d(j) = V(n-1,j);
		V(n-1,j) = 0.0;
	  }
	  V(n-1,n-1) = 1.0;
	  e(0) = 0.0;
	} 

	// Symmetric tridiagonal QL algorithm.
	inline void tql2(Matrix& V, Vector& d, Vector& e) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	  int n = d.size();
	  for (int i = 1; i < n; i++) {
		e(i-1) = e(i);
	  }
	  e(n-1) = 0.0;

	  double f = 0.0;
	  double tst1 = 0.0;
	  double eps = pow(2.0,-52.0);
	  for (int l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = EIG3_MAX(tst1,fabs(d(l)) + fabs(e(l)));
		int m = l;
		while (m < n) {
		  if (fabs(e(m)) <= eps*tst1) {
			break;
		  }
		  m++;
		}

		// If m == l, d(l) is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
		  int iter = 0;
		  do {
			iter = iter + 1;  // (Could check iteration count here.)

			// Compute implicit shift

			double g = d(l);
			double p = (d(l+1) - g) / (2.0 * e(l));
			double r = hypot2(p,1.0);
			if (p < 0) {
			  r = -r;
			}
			d(l) = e(l) / (p + r);
			d(l+1) = e(l) * (p + r);
			double dl1 = d(l+1);
			double h = g - d(l);
			for (int i = l+2; i < n; i++) {
			  d(i) -= h;
			}
			f = f + h;

			// Implicit QL transformation.

			p = d(m);
			double c = 1.0;
			double c2 = c;
			double c3 = c;
			double el1 = e(l+1);
			double s = 0.0;
			double s2 = 0.0;
			for (int i = m-1; i >= l; i--) {
			  c3 = c2;
			  c2 = c;
			  s2 = s;
			  g = c * e(i);
			  h = c * p;
			  r = hypot2(p,e(i));
			  e(i+1) = s * r;
			  s = e(i) / r;
			  c = p / r;
			  p = c * d(i) - s * g;
			  d(i+1) = h + s * (c * g + s * d(i));

			  // Accumulate transformation.

			  for (int k = 0; k < n; k++) {
				h = V(k,i+1);
				V(k,i+1) = s * V(k,i) + c * h;
				V(k,i) = c * V(k,i) - s * h;
			  }
			}
			p = -s * s2 * c3 * el1 * e(l) / dl1;
			e(l) = s * p;
			d(l) = c * p;

			// Check for convergence.

		  } while (fabs(e(l)) > eps*tst1);
		}
		d(l) = d(l) + f;
		e(l) = 0.0;
	  }
  
	  // Sort eigenvalues and corresponding vectors.

	  for (int i = 0; i < n-1; i++) {
		int k = i;
		double p = d(i);
		for (int j = i+1; j < n; j++) {
		  if (d(j) < p) {
			k = j;
			p = d(j);
		  }
		}
		if (k != i) {
		  d(k) = d(i);
		  d(i) = p;
		  for (int j = 0; j < n; j++) {
			p = V(j,i);
			V(j,i) = V(j,k);
			V(j,k) = p;
		  }
		}
	  }
	}

	inline void eigen_decomposition_sym_3x3(const Matrix& A, Matrix& V, Vector& d) {
	  int n = d.size();
	  Vector e(3);
	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  V(i,j) = A(i,j);
		}
	  }
	  tred2(V, d, e);
	  tql2(V, d, e);
	}

	// =====================================================================================================
	//
	// C.LAW
	//
	// =====================================================================================================

	DamageTC3DLaw::DamageTC3DLaw() 
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_rt(0.0)
		, m_rt_converged(0.0)
		, m_rc(0.0)
		, m_rc_converged(0.0)
		, m_damage_t(0.0)
		, m_damage_c(0.0)
		, m_lch(0.0)
		, m_lch_multiplier(1.0)
		, m_initial_strain()
#ifdef DAM_TC_3D_IMPLEX
		, m_rt_converged_old(0.0)
		, m_rc_converged_old(0.0)
		, m_strain()
		, m_strain_converged()
		, m_dTime_n(0.0)
		, m_dTime_n_converged(0.0)
		, m_rt_impl_temp(0.0)
		, m_rc_impl_temp(0.0)
#endif // DAM_TC_3D_IMPLEX
		, m_error_code(0.0)
	{
	}

	ConstitutiveLaw::Pointer DamageTC3DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer( new DamageTC3DLaw() );
	}

	DamageTC3DLaw::SizeType DamageTC3DLaw::WorkingSpaceDimension()
	{
		return 3;
	}

	DamageTC3DLaw::SizeType DamageTC3DLaw::GetStrainSize()
	{
		return 6;
	}

	bool DamageTC3DLaw::Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == DAMAGE_C)
			return true;
		if(rThisVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			return true;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			return true;
		return false;
	}

	bool DamageTC3DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_STRAIN)
			return true;
		return false;
	}

	bool DamageTC3DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool DamageTC3DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool DamageTC3DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	double& DamageTC3DLaw::GetValue(
		const Variable<double>& rThisVariable, 
		double& rValue)
	{
		rValue = 0.0;
		if(rThisVariable == DAMAGE_T)
			rValue = m_damage_t;
		else if(rThisVariable == DAMAGE_C)
			rValue = m_damage_c;
		else if(rThisVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			rValue = m_lch_multiplier;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		else if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		return rValue;
	}

	Vector& DamageTC3DLaw::GetValue(
		const Variable<Vector>& rThisVariable, 
		Vector& rValue)
	{
		if(rThisVariable == INITIAL_STRAIN) {
			if(rValue.size() != m_initial_strain.size())
				rValue.resize(m_initial_strain.size());
			noalias(rValue) = m_initial_strain;
		}
		return rValue;
	}

	Matrix& DamageTC3DLaw::GetValue(
		const Variable<Matrix>& rThisVariable, 
		Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & DamageTC3DLaw::GetValue(
		const Variable<array_1d<double, 3 > >& rVariable, 
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & DamageTC3DLaw::GetValue(
		const Variable<array_1d<double, 6 > >& rVariable, 
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void DamageTC3DLaw::SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable==TEMPERATURE)
			m_verbose = true;
		if(rVariable == DAMAGE_T)
			m_damage_t = rValue;
		else if(rVariable == DAMAGE_C)
			m_damage_c = rValue;
		else if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			m_lch_multiplier = rValue;
	}

	void DamageTC3DLaw::SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue, 
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == INITIAL_STRAIN) {
			if(rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
	}

	void DamageTC3DLaw::SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue, 
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTC3DLaw::SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTC3DLaw::SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool DamageTC3DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if( !rMaterialProperties.Has(YOUNG_MODULUS) ) return false;
		if( !rMaterialProperties.Has(POISSON_RATIO) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRESS_T_0) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_0) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_P) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_R) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRAIN_C_P) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_C) ) return false;
		if( !rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER) ) return false;
		if( !rMaterialProperties.Has(VISCOSITY) ) return false;
		return true;
	}

	DamageTC3DLaw::StrainMeasure DamageTC3DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	DamageTC3DLaw::StressMeasure DamageTC3DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool DamageTC3DLaw::IsIncremental()
	{
		return false;
	}

	void DamageTC3DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if(!m_initialized)
		{
			// random imperfections
			double impf = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(
				rElementGeometry,rShapeFunctionsValues);

			m_rt             = rMaterialProperties[DAMAGE_STRESS_T_0]*impf;
			m_rt_converged   = m_rt;
			m_rc             = rMaterialProperties[DAMAGE_STRESS_C_0]*impf;
			m_rc_converged   = m_rc;
			m_damage_t       = 0.0;
			m_damage_c       = 0.0;
#ifndef DAM_TC_OPTIMIZE_LCH
			m_lch            = rElementGeometry.Length();
#else
			m_lch            = get_optimized_lch(rElementGeometry);
#endif // !DAM_TC_OPTIMIZE_LCH
			m_lch_multiplier = 1.0;
			m_initial_strain = ZeroVector(this->GetStrainSize());
			m_initialized    = true;
#ifdef DAM_TC_3D_IMPLEX
			m_rt_converged_old = m_rt;
			m_rc_converged_old = m_rc;
			m_strain = ZeroVector(this->GetStrainSize());
			m_strain_converged = ZeroVector(this->GetStrainSize());
			m_dTime_n = 0.0;
			m_dTime_n_converged = 0.0;
#endif // DAM_TC_3D_IMPLEX

			// MOD
			/*double yavg = 0.0;
			for(unsigned int i=0; i<rElementGeometry.PointsNumber(); i++)
				yavg += rElementGeometry[i].Y0();
			yavg /= double(rElementGeometry.PointsNumber());
			if(yavg > 3600.0) {
				double fact = 2.0/3.0;
				m_rt            *= fact;
				m_rt_converged  *= fact;
				m_rc            *= fact;
				m_rc_converged  *= fact;
			}*/
			// END MOD

		}
		m_verbose = false;
		m_error_code = 0.0;
	}

	void DamageTC3DLaw::InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTC3DLaw::FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
#ifdef DAM_TC_3D_IMPLEX

		m_rt = m_rt_impl_temp;
		m_rc = m_rc_impl_temp;

		// move from n to n-1
		m_rt_converged_old  = m_rt_converged;
		m_rc_converged_old  = m_rc_converged;
		m_dTime_n_converged = m_dTime_n;
		m_strain_converged = m_strain;

#endif // DAM_TC_3D_IMPLEX

		// save converged values
		m_rt_converged = m_rt;
		m_rc_converged = m_rc;
	}

	void DamageTC3DLaw::InitializeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTC3DLaw::FinalizeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTC3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		const ProcessInfo&  pinfo = rValues.GetProcessInfo();
		const GeometryType& geom  = rValues.GetElementGeometry();
		const Properties&   props = rValues.GetMaterialProperties();

		const Vector& strain_vector       = rValues.GetStrainVector();
		Vector&       stress_vector       = rValues.GetStressVector();

		CalculationData data;

		bool verb = m_verbose;
		this->InitializeCalculationData(props, geom, rValues.GetShapeFunctionsValues(), pinfo, data);
		this->CalculateMaterialResponseInternal(strain_vector, stress_vector, data);
		m_verbose = false;

		/*if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
			if(props.Has(DAMAGE_SECANT_MATRIX)) {
				if(props[DAMAGE_SECANT_MATRIX] > 0) {
					size_t n = GetStrainSize();
					Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
					if(constitutive_matrix.size1() != n || constitutive_matrix.size2() != n)
						constitutive_matrix.resize(n, n);

					Matrix W( IdentityMatrix(3,3) );
					noalias(W) -= m_damage_t*data.PT;
					noalias(W) -= m_damage_c*data.PC;

					noalias(constitutive_matrix) = prod(W, data.C0);

					return;
				}
			}
		}*/

		//std::stringstream ss;

		if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR))
		{
			size_t n = GetStrainSize();

			// prepare constitutive matrix
			Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
			if(constitutive_matrix.size1() != n || constitutive_matrix.size2() != n)
				constitutive_matrix.resize(n, n);

			// save internal variables
			double save_rt = m_rt;
			double save_rc = m_rc;
			double save_dt = m_damage_t;
			double save_dc = m_damage_c;

			// perturbation parameter
			double h = 1.0E-9;

			// perturbed vectors
			Vector strain_bar(n);
			Vector S1(n);
			Vector S2(n);

			// apply perturbation to each strain component...
			for(size_t j = 0; j < n; j++)
			{
				h = std::max(1.0e-9, 1.0e-6*std::abs(strain_vector(j)));

				noalias(strain_bar) = strain_vector;

				/*strain_bar(j) = strain_vector(j) - h;
				this->CalculateMaterialResponseInternal(strain_bar, S1, data);

				strain_bar(j) = strain_vector(j) + h;
				this->CalculateMaterialResponseInternal(strain_bar, S2, data);

				for(size_t i = 0; i < n; i++)
					constitutive_matrix(i,j) = (S2(i) - S1(i))/(2.0*h);*/

				strain_bar(j) = strain_vector(j) + h;
				this->CalculateMaterialResponseInternal(strain_bar, S2, data);

				/*ss << MathHelpers::VectorToString(data.S, 4, std::scientific);
				ss << MathHelpers::VectorToString(data.Si, 4, std::scientific);
				ss << MathHelpers::VectorToString(data.ST, 4, std::scientific);
				ss << MathHelpers::VectorToString(data.SC, 4, std::scientific);
				ss << MathHelpers::VectorToString(S2, 4, std::scientific);
				ss << "****\n";*/

				for(size_t i = 0; i < n; i++)
					constitutive_matrix(i,j) = (S2(i) - stress_vector(i))/h;
			}

			// restore internal variables
			m_rt	   = save_rt;
			m_rc	   = save_rc;
			m_damage_t = save_dt;
			m_damage_c = save_dc;
		}

		/*ss << " ----------------\n";
		ss << "E:\n";
		ss << data.E << std::endl;
		ss << "STRAIN:\n";
		ss << MathHelpers::VectorToString(strain_vector, 4, std::scientific);
		ss << "STRESS:\n";
		ss << MathHelpers::VectorToString(stress_vector, 4, std::scientific);
		ss << "TANGENT:\n";
		ss << MathHelpers::MatrixToString(rValues.GetConstitutiveMatrix(), 4, std::scientific);
		ss << "C0:\n";
		ss << MathHelpers::MatrixToString(data.C0, 4, std::scientific);
		if(data.E < 50000.0) {
			std::cout << ss.str();
			exit(-1);
		}

		m_verbose = verb;*/
	}

	void DamageTC3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTC3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{

	}

	void DamageTC3DLaw::ResetMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		m_rt = 0.0;
		m_rt_converged = 0.0;
		m_rc = 0.0;
		m_rc_converged = 0.0;
		m_damage_t = 0.0;
		m_damage_c = 0.0;
		m_lch = 0.0;
		m_lch_multiplier = 1.0;
		m_initialized = false;
	}

	void DamageTC3DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int DamageTC3DLaw::Check(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if( !rMaterialProperties.Has(YOUNG_MODULUS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

		if( !rMaterialProperties.Has(POISSON_RATIO) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

		if( !rMaterialProperties.Has(DAMAGE_STRESS_T_0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRESS_T_0", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_T", "");

		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRESS_C_0", "");

		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_P) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRESS_C_P", "");

		if( !rMaterialProperties.Has(DAMAGE_STRESS_C_R) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRESS_C_R", "");

		if( !rMaterialProperties.Has(DAMAGE_STRAIN_C_P) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRAIN_C_P", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_C) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_C", "");
		
		if( !rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: BIAXIAL_COMPRESSION_MULTIPLIER", "");

		return 0;

		KRATOS_CATCH("");
	}

	void DamageTC3DLaw::CalculateMaterialResponse(const Vector& StrainVector,
			const Matrix& DeformationGradient,
			Vector& StressVector,
			Matrix& AlgorithmicTangent,
			const ProcessInfo& rCurrentProcessInfo,
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			bool CalculateStresses,
			int CalculateTangent,
			bool SaveInternalVariables)
	{
		ConstitutiveLaw::Parameters parameters(rElementGeometry, rMaterialProperties, rCurrentProcessInfo);
		Vector E(StrainVector);
		parameters.SetStrainVector( E );
		parameters.SetStressVector( StressVector );
		parameters.SetConstitutiveMatrix( AlgorithmicTangent );
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateStresses);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateTangent);
		double detF = 1.0;
		Matrix F(IdentityMatrix(3,3));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);
		parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
		this->CalculateMaterialResponseCauchy(parameters);
	}

	void DamageTC3DLaw::InitializeCalculationData(const Properties& props, 
												  const GeometryType& geom,
												  const Vector& N,
												  const ProcessInfo& pinfo,
												  CalculationData& data)
	{
		// random imperfections
		double impf = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(geom, N);

		// elasticity
		data.E   = props[YOUNG_MODULUS]*impf;
		data.nu  = props[POISSON_RATIO];
		this->CalculateElasticityMatrix(data);

		// tension
		data.ft = props[DAMAGE_STRESS_T_0]*impf;
		data.Gt = props[FRACTURE_ENERGY_T];
		data.m2  = props.Has(DAMAGE_TENSILE_SURFACE_S1) ? props[DAMAGE_TENSILE_SURFACE_S1] : 1.0;
		data.m2  = std::min(std::max(data.m2,0.0),1.0);

		// compression
		data.fc0 = props[DAMAGE_STRESS_C_0]*impf;
		data.fcp = props[DAMAGE_STRESS_C_P]*impf;
		data.fcr = props[DAMAGE_STRESS_C_R];
		data.ep  = props[DAMAGE_STRAIN_C_P]*impf;
		data.c1  = props.Has(DAMAGE_COMPRESSIVE_LAW_C1) ? props[DAMAGE_COMPRESSIVE_LAW_C1] : 0.65;
		data.c2  = props.Has(DAMAGE_COMPRESSIVE_LAW_C2) ? props[DAMAGE_COMPRESSIVE_LAW_C2] : 0.50;
		data.c3  = props.Has(DAMAGE_COMPRESSIVE_LAW_C3) ? props[DAMAGE_COMPRESSIVE_LAW_C3] : 1.50;
		data.Gc  = props[FRACTURE_ENERGY_C];
		data.bm  = props[BIAXIAL_COMPRESSION_MULTIPLIER];
		data.m1  = props.Has(SHEAR_COMPRESSION_REDUCTION) ? props[SHEAR_COMPRESSION_REDUCTION] : 0.5;
		data.m1  = std::min(std::max(data.m1,0.0),1.0);

		// MOD
		//double yavg = 0.0;
		//for(unsigned int i=0; i<geom.PointsNumber(); i++)
		//	yavg += geom[i].Y0();
		//yavg /= double(geom.PointsNumber());
		//if(yavg > 3600.0) {
		//	double fact = 2.0/3.0;
		//	//data.E*=fact;
		//	data.fc0*=fact;
		//	data.fcp*=fact;
		//	data.fcm*=fact;
		//	data.fcr*=fact;
		//	data.Gc*=fact;
		//	data.ft*=fact;
		//	data.Gt*=fact;
		//}
		// END MOD

		// lubliner parameter for compressive states
		data.Kc = 2.0/3.0;
		if(props.Has(LUBLINER_SURFACE_PARAM_KC)) {
			data.Kc = props[LUBLINER_SURFACE_PARAM_KC];
		}
		
		// effective stress data
		data.S.resize(6,false);
		data.Si.resize(3,false);
		data.ST.resize(6,false);
		data.SC.resize(6,false);
		//data.PT.resize(6,6,false);
		//data.PC.resize(6,6,false);

		// misc
		data.eta = props[VISCOSITY];
		data.lch = m_lch * m_lch_multiplier;

		data.dTime = pinfo[DELTA_TIME];
		if(data.dTime > 0.0 && data.eta > 0.0)
		{
			data.rate_coeff_1 = data.eta/(data.eta+data.dTime);
			data.rate_coeff_2 = data.dTime/(data.eta+data.dTime);
		}
		else
		{
			data.rate_coeff_1 = 0.0;
			data.rate_coeff_2 = 1.0;
		}

		// tensile model:
		/*
		0 = lubliner (default)
		1 = rankine 
		*/
		data.tensile_damage_model = props.Has(DAMAGE_TENSILE_MODEL) ? props[DAMAGE_TENSILE_MODEL] : 0;
	}

	void DamageTC3DLaw::CalculateElasticityMatrix(CalculationData& data)
	{
		if(data.C0.size1() != 6 || data.C0.size2() != 6)
			data.C0.resize(6,6,false);
		data.C0.clear();
		
		double c1 = data.E / (1.0 - data.nu * data.nu);
		double c2 = c1 * data.nu;
		double c3 = c1 * (1.0 - data.nu) / 2.0;

		data.C0(0, 0) = (data.E*(1.0 - data.nu) / ((1.0 + data.nu)*(1.0 - 2.0*data.nu)));
		data.C0(1, 1) = data.C0(0, 0);
		data.C0(2, 2) = data.C0(0, 0);

		data.C0(3, 3) = data.C0(0, 0)*(1.0 - 2.0 * data.nu) / (2.0*(1.0 - data.nu));
		data.C0(4, 4) = data.C0(3, 3);
		data.C0(5, 5) = data.C0(3, 3);

		data.C0(0, 1) = data.C0(0, 0)*data.nu / (1.0 - data.nu);
		data.C0(1, 0) = data.C0(0, 1);

		data.C0(0, 2) = data.C0(0, 1);
		data.C0(2, 0) = data.C0(0, 1);

		data.C0(1, 2) = data.C0(0, 1);
		data.C0(2, 1) = data.C0(0, 1);
	}

	void DamageTC3DLaw::TensionCompressionSplit(CalculationData& data)
	{
		const Vector& S  = data.S;
		Vector&       Si = data.Si;
		Vector&       ST = data.ST;
		Vector&       SC = data.SC;

		if(Si.size() != 3) Si.resize(3,false);
		if(ST.size() != 6) ST.resize(6,false);
		if(SC.size() != 6) SC.resize(6,false);

		ST.clear();
		SC.clear();

		array_1d<double,3> v0;
		array_1d<double,3> v1;
		array_1d<double,3> v2;
		array_1d<double,3> vtemp;
		
		double p1 = S(3)*S(3) + S(4)*S(4) + S(5)*S(5);
		/*if(p1 <= DAM_TC_PREC)*/
		if(p1 == 0.0)
		{
			v0(0)=1.0; v0(1)=0.0; v0(2)=0.0;
			v1(0)=0.0; v1(1)=1.0; v1(2)=0.0;
			v2(0)=0.0; v2(1)=0.0; v2(2)=1.0;
			Si(0)=S(0);
			Si(1)=S(1);
			Si(2)=S(2);
		}
		else
		{
			Matrix V(3,3);
			Matrix A(3,3);
			A(0,0) = S(0);
			A(1,1) = S(1);
			A(2,2) = S(2);
			A(0,1) = S(3);
			A(1,2) = S(4);
			A(0,2) = S(5);
			A(1,0) = S(3);
			A(2,1) = S(4);
			A(2,0) = S(5);

			eigen_decomposition_sym_3x3(A,V,Si);
			for(int i=0; i<3; i++) {
				v0(i) = V(i,0);
				v1(i) = V(i,1);
				v2(i) = V(i,2);
				// watch out! the algo above assumes eigenvectors on rows!!! bloodclat!!!!
			}
			v0 /= norm_2(v0);
			v1 /= norm_2(v1);
			v2 /= norm_2(v2);
		}
		if (Si(0) < Si(1)) {
			dam_tc_swap(Si(0),Si(1));
			dam_tc_swap(v0,v1,vtemp);
		}
		if (Si(1) < Si(2)) {
			dam_tc_swap(Si(1),Si(2));
			dam_tc_swap(v1,v2,vtemp);
		}
		if (Si(0) < Si(1)) {
			dam_tc_swap(Si(0),Si(1));
			dam_tc_swap(v0,v1,vtemp);
		}
		if(Si(0)>0.0)
			dam_tc_add_dyad(v0, Si(0), ST);
		else if(Si(0)<0.0)
			dam_tc_add_dyad(v0, Si(0), SC);
		if(Si(1)>0.0)
			dam_tc_add_dyad(v1, Si(1), ST);
		else if(Si(1)<0.0)
			dam_tc_add_dyad(v1, Si(1), SC);
		if(Si(2)>0.0)
			dam_tc_add_dyad(v2, Si(2), ST);
		else if(Si(2)<0.0)
			dam_tc_add_dyad(v2, Si(2), SC);
	}
	
	void DamageTC3DLaw::TrialEquivalentStressTension(CalculationData& data, double& rt_trial)
	{
		if(data.tensile_damage_model == 0) // LUBLINER
		{
			double fc = data.fcp*data.m2;
			double fb = fc * data.bm;
			double alpha = (fb-fc)/(2.0*fb-fc);
			double ft = data.ft;
			double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);
			double I1 = data.Si(0) + data.Si(1) + data.Si(2);
			array_1d<double,3>Sdev;
			for(size_t i=0; i < 3; i++)
				Sdev(i) = data.Si(i)-I1/3.0;
			double J2 = 0.5*inner_prod(Sdev,Sdev);
			double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
			double smax_alg = data.Si(0);
			double smax = std::max( smax_alg,0.0);
			double smin = std::max(-smax_alg,0.0);
			rt_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + beta*smax - gamma*smin) /fc*ft;
		}
		else if(data.tensile_damage_model == 1) // RANKINE
		{
			rt_trial = std::max(std::max(data.Si(0), data.Si(1)), 0.0);
		}
		else // NO DAMAGE IN TENSION
		{
			rt_trial = 0.0;
		}

		double sign_check = data.Si(0) > 0.0 ? 1.0 : 0.0;
		rt_trial *= sign_check;

		//if( (DAM_TC_SIGN(data.Si(0)) + DAM_TC_SIGN(data.Si(1)) + DAM_TC_SIGN(data.Si(2))) < 0.0 )
		//{
		//	rt_trial = 0.0;
		//}
		//else
		//{

		//	if(data.Si(0) < 0.0) {
		//		double daux = 1.0E-7*std::abs(data.Si(0));
		//		if(data.Si(1) < daux && data.Si(2) < daux) {
		//			rt_trial = 0.0;
		//			return;
		//		}
		//	}
		//	if(data.Si(1) < 0.0) {
		//		double daux = 1.0E-7*std::abs(data.Si(1));
		//		if(data.Si(0) < daux && data.Si(2) < daux) {
		//			rt_trial = 0.0;
		//			return;
		//		}
		//	}
		//	if(data.Si(2) < 0.0) {
		//		double daux = 1.0E-7*std::abs(data.Si(2));
		//		if(data.Si(0) < daux && data.Si(1) < daux) {
		//			rt_trial = 0.0;
		//			return;
		//		}
		//	}


		//	if(data.tensile_damage_model == 0) // LUBLINER
		//	{
		//		double fc = data.fcp*data.m2;
		//		double fb = fc * data.bm;
		//		double alpha = (fb-fc)/(2.0*fb-fc);
		//		double ft = data.ft;
		//		double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);
		//		double I1 = data.Si(0) + data.Si(1) + data.Si(2);
		//		array_1d<double,3>Sdev;
		//		for(size_t i=0; i < 3; i++)
		//			Sdev(i) = data.Si(i)-I1/3.0;
		//		double J2 = 0.5*inner_prod(Sdev,Sdev);
		//		double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
		//		double smax_alg = data.Si(0);
		//		double smax = std::max( smax_alg,0.0);
		//		double smin = std::max(-smax_alg,0.0);
		//		rt_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + beta*smax - gamma*smin) /fc*ft;
		//	}
		//	else if(data.tensile_damage_model == 1) // RANKINE
		//	{
		//		rt_trial = std::max(std::max(data.Si(0), data.Si(1)), 0.0);
		//	}
		//	else // NO DAMAGE IN TENSION
		//	{
		//		rt_trial = 0.0;
		//	}

		//	//if(data.tensile_damage_model == 1) { // rankine
		//	//	rt_trial = std::max(data.Si(0),0.0);
		//	//	return;
		//	//}

		//	//double fc = data.fcp*data.m2;
		//	//double fb = fc * data.bm;
		//	//double alpha = (fb-fc)/(2.0*fb-fc);
		//	//double ft = data.ft;

		//	//double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);
		//	//
		//	//double I1 = data.Si(0) + data.Si(1) + data.Si(2);
		//	//array_1d<double,3>Sdev;
		//	//for(size_t i=0; i < 3; i++)
		//	//	Sdev(i) = data.Si(i)-I1/3.0;
		//	//double J2 = 0.5*inner_prod(Sdev,Sdev);

		//	//double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
		//	//double smax_alg = data.Si(0);
		//	//double smax = std::max( smax_alg,0.0);
		//	//double smin = std::max(-smax_alg,0.0);

		//	//rt_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + beta*smax - gamma*smin) /fc*ft;
		//}
	}

	void DamageTC3DLaw::TrialEquivalentStressCompression(CalculationData& data, double& rc_trial)
	{
		double fc = data.fc0;
		double fb = fc * data.bm;
		double alpha = (fb-fc)/(2.0*fb-fc);
		double ft = data.ft*data.fc0/data.fcp;

		double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);

		double I1 = data.Si(0) + data.Si(1) + data.Si(2);
		array_1d<double,3>Sdev;
		for(size_t i=0; i < 3; i++)
			Sdev(i) = data.Si(i)-I1/3.0;
		double J2 = 0.5*inner_prod(Sdev,Sdev);

		double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
		double smax_alg = data.Si(0);
		double smax = std::max( smax_alg,0.0);
		double smin = std::max(-smax_alg,0.0);
		rc_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + data.m1*beta*smax - gamma*smin);

		double sign_check = data.Si(2) < 0.0 ? 1.0 : 0.0;
		rc_trial *= sign_check;

		/*if( (DAM_TC_SIGN(data.Si(0)) + DAM_TC_SIGN(data.Si(1)) + DAM_TC_SIGN(data.Si(2))) > 0.0 )
		{
			rc_trial = 0.0;
		}
		else
		{
			if(data.Si(0) > 0.0) {
				double daux = 1.0E-7*std::abs(data.Si(0));
				if(data.Si(1) < daux && data.Si(2) < daux) {
					rc_trial = 0.0;
					return;
				}
			}
			if(data.Si(1) > 0.0) {
				double daux = 1.0E-7*std::abs(data.Si(1));
				if(data.Si(0) < daux && data.Si(2) < daux) {
					rc_trial = 0.0;
					return;
				}
			}
			if(data.Si(2) > 0.0) {
				double daux = 1.0E-7*std::abs(data.Si(2));
				if(data.Si(0) < daux && data.Si(1) < daux) {
					rc_trial = 0.0;
					return;
				}
			}
			
			double fc = data.fc0;
			double fb = fc * data.bm;
			double alpha = (fb-fc)/(2.0*fb-fc);
			double ft = data.ft*data.fc0/data.fcp;

			double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);

			double I1 = data.Si(0) + data.Si(1) + data.Si(2);
			array_1d<double,3>Sdev;
			for(size_t i=0; i < 3; i++)
				Sdev(i) = data.Si(i)-I1/3.0;
			double J2 = 0.5*inner_prod(Sdev,Sdev);

			double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
			double smax_alg = data.Si(0);
			double smax = std::max( smax_alg,0.0);
			double smin = std::max(-smax_alg,0.0);
			rc_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + data.m1*beta*smax - gamma*smin);
		}*/
	}

	void DamageTC3DLaw::CalculateDamageTension(CalculationData& data, double rt, double& dt)
	{
		if(rt <= data.ft)
		{
			dt = 0.0;
		}
		else
		{
			double lch = data.lch;
			double E   = data.E;
			double ft  = data.ft;
			double G   = data.Gt;
			double r0  = ft;
			double lt  = 2.0*E*G/ft/ft;

			if(lch >= lt) 
			{
				std::stringstream ss;
				ss << "FRACTURE_ENERGY_T is to low:  2*E*Gt/(ft*ft) = " << lt 
					<< ",   Characteristic Length = " << lch << std::endl;
				std::cout << ss.str();
				exit(-1);
			}

			double Hs  = lch/(lt-lch);
			double A   = 2.0*Hs; 
			dt         = 1.0-r0/rt*std::exp(A*(1.0-rt/r0));

			double rmin=1.0e-2*ft;
			if((1.0-dt)*rt<rmin)
				dt = 1.0-rmin/rt;
		}
	}

	void DamageTC3DLaw::CalculateDamageCompression(CalculationData& data, double rc, double& dc)
	{
		if(rc <= data.fc0)
		{
			dc = 0.0;
		}
		else
		{
			// extract material parameters
			double E  = data.E;
			double s0 = data.fc0;
			double sp = data.fcp;
			double sr = data.fcr;
			double ep = data.ep;
			double c1 = data.c1;
			double c2 = data.c2;
			double c3 = data.c3;
			double Gc = data.Gc / data.lch;

			// auto-computation of other parameters
			double sk    = sr + (sp-sr)*c1;
			double e0    = s0/E;
			double alpha = 2.0*(ep-sp/E);
			double ej    = ep + alpha*c2;
			double ek    = ej + alpha*(1-c2);
			double er    = (ek-ej)/(sp-sk)*(sp-sr)+ej;
			double eu    = er*c3;

			// regularization
			double G_bar   = b3_calc_G( sp,sk,sr,ep,ej,ek,er,eu );
			double G1      = sp*ep/2.0;
			double stretch = (Gc-G1)/(G_bar-G1)-1.0;
			if(stretch <= -1.0) 
			{
				std::stringstream ss;
				ss << "Damage TC Error: Compressive fracture energy is too low" << std::endl;
				ss << "Input Gc/lch = " << Gc << std::endl;
				ss << "Minimum Gc to avoid constitutive snap-back = " << G1 << std::endl;
				std::cout << ss.str();
				exit(-1);
			}
			b3_stretch( stretch,ep,ej,ek,er,eu );

			// current abscissa
			double xi = rc/E;

			// compute damage
			double s = rc;
			if(xi <= ep)
				s = b3_eval_bezier(xi,e0,sp/E,ep,s0,sp,sp);
			else if(xi <= ek)
				s = b3_eval_bezier(xi,ep,ej,ek,sp,sp,sk);
			else if(xi <= eu)
				s = b3_eval_bezier(xi,ek,er,eu,sk,sr,sr);
			else
				s = sr;
			dc = 1.0-s/rc;
		}
	}

	void DamageTC3DLaw::CalculateMaterialResponseInternal(const Vector& strain_vector,
																	 Vector& stress_vector,
																	 CalculationData& data)
	{
		size_t strain_size = this->GetStrainSize();

		if(stress_vector.size() != strain_size)
			stress_vector.resize(strain_size,false);

		// set up coefficients for the rate-dependent model

		double rate_coeff_1 = data.rate_coeff_1;
		double rate_coeff_2 = data.rate_coeff_2;

		// elastic predictor + stress decomposition (tension-compression split)

		m_rt = m_rt_converged;
		m_rc = m_rc_converged;

		noalias(data.S) = prod(data.C0, strain_vector - m_initial_strain);

		for(size_t i = 0; i < 6; i++)
			if(std::abs(data.S(i)) < DAM_TC_PREC) data.S(i) = 0.0;

		this->TensionCompressionSplit(data);

		// compute the equivalent stress measures

		double rt_trial;
		double rc_trial;
		this->TrialEquivalentStressTension(data, rt_trial);
		this->TrialEquivalentStressCompression(data, rc_trial);

		// damage update

#ifdef DAM_TC_3D_IMPLEX

		// time factor
		double time_factor = 0.0;
		if(m_dTime_n_converged>0.0) time_factor = data.dTime/m_dTime_n_converged;
		m_dTime_n = data.dTime;

		if(true)//data.dTime < 1.0E-4)
		{
			// IMPLEX INTEGRATION

			// explicit evaluation
			m_rt = m_rt_converged + time_factor * (m_rt_converged-m_rt_converged_old);
			m_rc = m_rc_converged + time_factor * (m_rc_converged-m_rc_converged_old);

			// save implicit variables for the finalize_solution_step
			double rt_impl = m_rt_converged;
			double rc_impl = m_rc_converged;
			if(rt_trial > rt_impl)
				rt_impl = rate_coeff_1*rt_impl + rate_coeff_2*rt_trial;
			if(rc_trial > rc_impl)
				rc_impl = rate_coeff_1*rc_impl + rate_coeff_2*rc_trial;
			m_rt_impl_temp = rt_impl;
			m_rc_impl_temp = rc_impl;

			// new damage variables (explicit)

			this->CalculateDamageTension(data, m_rt, m_damage_t);
			this->CalculateDamageCompression(data, m_rc, m_damage_c);
		}
		else
		{
			// IMPLICIT INTEGRATION
			if(rt_trial > m_rt)
				m_rt = rate_coeff_1*m_rt + rate_coeff_2*rt_trial;
			this->CalculateDamageTension(data, m_rt, m_damage_t);

			if(rc_trial > m_rc)
				m_rc = rate_coeff_1*m_rc + rate_coeff_2*rc_trial;
			this->CalculateDamageCompression(data, m_rc, m_damage_c);

			m_rt_impl_temp = m_rt;
			m_rc_impl_temp = m_rc;
		}

		// check explicit error
		m_error_code = 0.0;

		/*double dam_t_old, dam_c_old;
		this->CalculateDamageTension(data, m_rt_converged, dam_t_old);
		this->CalculateDamageCompression(data, m_rc_converged, dam_c_old);
		double d1 = m_damage_t - dam_t_old;
		double d2 = m_damage_c - dam_c_old;
		double dd = std::max(d1,d2);
		if(dd > 0.01) {
			m_error_code = -1.0;
		}*/
#else

		if(rt_trial > m_rt)
			m_rt = rate_coeff_1*m_rt + rate_coeff_2*rt_trial;
		this->CalculateDamageTension(data, m_rt, m_damage_t);

		if(rc_trial > m_rc)
			m_rc = rate_coeff_1*m_rc + rate_coeff_2*rc_trial;
		this->CalculateDamageCompression(data, m_rc, m_damage_c);

#endif // DAM_TC_3D_IMPLEX

		// calculation of stress tensor
		noalias(stress_vector)  = (1.0 - m_damage_t)*data.ST;
		noalias(stress_vector) += (1.0 - m_damage_c)*data.SC;
	}


} /* namespace Kratos.*/
