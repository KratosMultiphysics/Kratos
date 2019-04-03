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

#include "damage_tc_plane_stress_2d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"
#include "custom_utilities/imperfection_utilities.h"
#include "includes/variables.h"

#define DAM_TC_PREC 0.0

#define DAM_TC_SIGN(X) (X == 0.0 ? 0.0 : ( X > 0.0 ? 1.0 : -1.0 ))
#define DAM_TC_HVS(X)  ( X > 0.0 ? 1.0 : ( X < 0.0 ? -1.0 : 0.5) )
#define DAM_TC_POS(X)  ( X > 0.0 ? X : 0.0 )

#define DAM_TC_OPTIMIZE_LCH

namespace Kratos
{

	// todo: remove this. only used for very deformed rectangular geometries

#ifdef DAM_TC_OPTIMIZE_LCH

	inline double get_optimized_lch(const Element::GeometryType& geom)
	{
		double lch = geom.Length();
		if(geom.WorkingSpaceDimension() == 2) {
			if(geom.PointsNumber() == 4) {
				// 2D 4N
				double ax = (geom[0].X0()+geom[3].X0())/2.0;
				double ay = (geom[0].Y0()+geom[3].Y0())/2.0;
				double bx = (geom[1].X0()+geom[2].X0())/2.0;
				double by = (geom[1].Y0()+geom[2].Y0())/2.0;
				double cx = (geom[0].X0()+geom[1].X0())/2.0;
				double cy = (geom[0].Y0()+geom[1].Y0())/2.0;
				double dx = (geom[2].X0()+geom[3].X0())/2.0;
				double dy = (geom[2].Y0()+geom[3].Y0())/2.0;
				double v1x = bx-ax;
				double v1y = by-ay;
				double v2x = dx-cx;
				double v2y = dy-cy;
				double lx = std::sqrt(v1x*v1x+v1y*v1y);
				double ly = std::sqrt(v2x*v2x+v2y*v2y);
				lch = std::min(lx,ly);
			}
		}
		else if(geom.WorkingSpaceDimension() == 3) {
			if(geom.PointsNumber() == 4) {
				// 3D 4N
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

	// TODO: move this into a BIFURCATION_ANALYSIS_UTILS_H

	inline double det_acoustic_tensor_2d(const Matrix& C, double a)
	{
		double sinA = std::sin(a);
		double cosA = std::cos(a);
		double a1 = sinA*sinA;
		double a2 = cosA*cosA;
		double a3 = C(1, 2);
		double a4 = C(2, 0);
		double a5 = C(2, 1);
		double a6 = C(0, 2);
		double a7 = C(1, 0);
		double a8 = C(0, 1);
		double a9 = C(1, 1);
		double a10 = C(0, 0);
		double a11 = C(2, 2);
		double a12 = a2*a2;
		double a13 = a1*a1;
		double detQ =
				(a1*a4*a8*cosA + a2*a3*a8*cosA - a1*a3*a10*cosA +
				a1*a6*a7*cosA + a2*a5*a7*cosA - a2*a4*a9*cosA -
				a1*a5*a10*cosA - a2*a6*a9*cosA)*sinA +
				a9*a11*a12 - a4*a6*a13 - a3*a5*a12 + a10*a11*a13 +
				a1*a2*a3*a4 + a1*a2*a5*a6 - a1*a2*a7*a8 -
				a1*a2*a7*a11 - a1*a2*a8*a11 + a1*a2*a9*a10;
		return detQ;
	}

	inline double min_det_acoustic_tensor_2d(const Matrix& C, double& angle)
	{
		double dq = std::numeric_limits<double>::max();
		angle = 0.0;
		for(int i = 0; i <= 180.0; i++)
		{
			double a = double(i)/180.0*Globals::Pi;
			double dq_trial = det_acoustic_tensor_2d(C,a);
			if(dq_trial < dq) {
				dq = dq_trial;
				angle = a;
			}
		}
		return dq;
	}

	// TODO: move this into a B3_HARDENING_LAW_H

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

	inline void bz_int(double x0,double x1,double x2,
		               double y0,double y1,double y2,
					   double Sx,double Sy,
					   double& xb, double& yb)
	{
		double k = -Sy/Sx;
		double a = k*(x0-2.0*x1+x2)+y0-2.0*y1+y2;

		if(std::abs(a)<1.0E-12) {
			x1 = x1+1.0E-6*(x2-x0);
			a = k*(x0-2.0*x1+x2)+y0-2.0*y1+y2;
		}

		double b = -2.0*(k*(x0-x1)+y0-y1);
		double c = k*(x0)+y0;

		double delta=b*b-4.0*a*c;
		double t1 = (-b+std::sqrt(delta))/(2.0*a);
		double t2 = (-b-std::sqrt(delta))/(2.0*a);

		double tol = 1.0e-5;
		double t=0.0;
		if(t1>=-tol && t1<=1.0+tol) {
			t = t1;
		}
		else {
			if(t2>=-tol && t2<=1.0+tol) {
				t = t2;
			}
			else {
				xb=0.0;
				yb=0.0;
				return;
			}
		}
		xb = (x0-2.0*x1+x2)*t*t+2.0*(x1-x0)*t+x0;
		yb = (y0-2.0*y1+y2)*t*t+2.0*(y1-y0)*t+y0;
	}

	inline double bz_2_f(double Sx,double Sy,
		                 double x0,double x1,double x2,
						 double y0,double y1,double y2,double ft)
	{
		double xb,yb;
		bz_int(x0,x1,x2,y0,y1,y2, Sx,Sy, xb,yb);
		double L0 = std::sqrt(xb*xb+yb*yb);
		double L1 = std::sqrt(Sx*Sx+Sy*Sy);
		double F = 0.0;
		if(L0>0.0)
			F=L1/L0*ft;//-ft;
		return F;
	}

	// todo: move this into a DAMAGE_INCREMENTAL_FRACTURE_ENERGY_REGULARIZATION_H

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

#define REG2_DIS_VEC_T std::vector<double>
#define REG2_DIS_VEC_C std::vector<double>

	template<class T>
	inline unsigned int reg2_find_stop_id(const T& x, const T& y, double xstop)
	{
		unsigned int id = -1;
		unsigned int n = x.size();
		for(unsigned int i=0; i<n-1; i++) {
			if(y[i+1]<y[i]) {
				id = i;
				break;
			}
		}
		if(id == -1) {
			id = n-1;
		}
		else {
			for(unsigned int i=0; i<n; i++) {
				if(x[i]>=xstop) {
					id = i;
					break;
				}
			}
		}
		return id;
	}

	inline void reg2_discrete_t_law(double E, double G, double s, REG2_DIS_VEC_T& x, REG2_DIS_VEC_T& y)
	{
		double h = 1.0;
		double lt=2.0*E*G/s/s;
		double Hs = h/(lt-h);
		double A  = 2.0*Hs;
		double r0 = s;

		double smin = 1.0e-3*s;
		double ru = -r0*(log(smin/r0)/A - 1.0);
		double xu = ru/E;
		double x0 = s/E;

		unsigned int n = 200;
		x.resize(n);
		y.resize(n);

		double x_incr = (xu-x0)/double(n - 1);
		for(unsigned int i=0; i<n; i++) {
			x[i] = x0 + double(i)*x_incr;
			double rt = x[i]*E;
			double d = 1.0-r0/rt*std::exp(A*(1.0-rt/r0));
			y[i] = (1-d)*rt;
		}
	}

	inline void reg2_discrete_c_law(const DamageTCPlaneStress2DLaw::CalculationData& data, double Gc,
									REG2_DIS_VEC_C& x, REG2_DIS_VEC_C& y)
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
		// auto-computation of other parameters
		double sk    = sr + (sp-sr)*c1;
		double e0    = s0/E;
		double alpha = 2.0*(ep-sp/E);
		double ej    = ep + alpha*c2;
		double ek    = ej + alpha*(1-c2);
		double er    = (ek-ej)/(sp-sk)*(sp-sr)+ej;
		double eu    = er*c3;

		// make the underneath area = Gc
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

		// generate lin-spaced abscissas
		unsigned int n1 = 41;
		unsigned int n2 = 80;
		unsigned int n3 = 80;
		unsigned int n  = n1+n2+n3;
		x.resize(n);
		y.resize(n);
		double dx1 = (ep-e0)/double(n1-1);
		double dx2 = (ek-ep)/double(n2);
		double dx3 = (eu-ek)/double(n3);
		unsigned int w=0;
		for(unsigned int i=0; i<n1; i++) x[w++]=e0+double(i  )*dx1;
		for(unsigned int i=0; i<n2; i++) x[w++]=ep+double(i+1)*dx2;
		for(unsigned int i=0; i<n3; i++) x[w++]=ek+double(i+1)*dx3;

		// generate curve
		for(unsigned int i=0; i<n; i++) {
			double xi = x[i];
			if(xi <= ep)
				y[i] = b3_eval_bezier(xi,e0,sp/E,ep,s0,sp,sp);
			else if(xi <= ek)
				y[i] = b3_eval_bezier(xi,ep,ej,ek,sp,sp,sk);
			else if(xi <= eu)
				y[i] = b3_eval_bezier(xi,ek,er,eu,sk,sr,sr);
			else
				y[i] = sr;
		}
	}

	template<class T>
	inline double reg2_interp_dis_curve(const T& x, const T& y, double ix)
	{
		double iy=0.0;
		unsigned int n = x.size();
		if(ix <= x[0]) {
			iy = y[0];
		}
		else {
			if(ix >= x[n-1]) {
				iy = y[n-1];
			}
			else {
				for(unsigned int i=0; i<n-1; i++) {
					double x1 = x[i];
					if(ix >= x1) {
						double x2 = x[i+1];
						if(ix < x2) {
							double y1 = y[i];
							double y2 = y[i+1];
							iy = (ix-x1)/(x2-x1)*(y2-y1)+y1;
						}
					}
				}
			}
		}
		return iy;
	}

	template<class T>
	inline void reg2_correct_snap_back(T& x)
	{
		unsigned int n = x.size();
		bool found = false;
		double last_incr=0.0;
		for(unsigned int i=1; i<n; i++) {
			if(x[i]<x[i-1]) {
				x[i] = x[i-1]+last_incr*1.0E-2;
				found=true;
			}
			else {
				if(!found)
					last_incr = x[i]-x[i-1];
			}
		}
	}

	template<class T>
	inline void reg2_regularize_from_point(T& x, const T& y, unsigned int pn, double lch_old, double lch_new)
	{
		//double xp = x[pn];
		//double yp = y[pn];
		//double Ed = yp/xp; // damaged young modulus
		//unsigned int n = x.size();
		//for(unsigned int i=pn; i<n; i++) {
		//	double xtot = x[i]; // total strain
		//	double xinl = xtot-y[i]/Ed; // inelastic strain
		//	double xinl_reg = xinl*lch_old/lch_new; // m_lch/lch = m_lch/(m_lch*m_lch_mult) = 1/m_lch_mult
		//	x[i] = xinl_reg + y[i]/Ed; // to regularized total strain
		//}

		double xp = x[pn];
		double yp = y[pn];
		unsigned int n = x.size();
		double G1 = xp*yp/2.0;
		double G2 = 0.0;
		for(unsigned int i = pn; i < n; i++)
			G2 += (x[i] - x[i-1])*(y[i]+y[i-1])/2.0;
		double G = G1+G2;
		double Gr = G/lch_new*lch_old;
		double S = (Gr-G1)/(G-G1)-1.0;
		if(S < -1.0)
			S = -0.999;
		for(unsigned int i = pn; i < n; i++)
			x[i] += (x[i]-xp)*S;
	}

#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

	// todo: move this into a EIGENVAL_UTILITIES_H

	/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public domain Java Matrix library JAMA. */

	#ifdef EIG3_MAX
	#undef EIG3_MAX
	#endif
	#define EIG3_MAX(a, b) ((a)>(b)?(a):(b))

	inline double hypot2(double x, double y) {
		return std::sqrt(x*x+y*y);
	}

	// Symmetric Householder reduction to tridiagonal form.
	template<class TMatrix, class TVector>
	inline void tred2(TMatrix& V, TVector& d, TVector& e) {

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
	template<class TMatrix, class TVector>
	inline void tql2(TMatrix& V, TVector& d, TVector& e) {

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

	inline void eigen_decomposition_sym_2x2(const Matrix& A, Matrix& V, Vector& d) {
	  int n = d.size();
	  Vector e(2);
	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  V(i,j) = A(i,j);
		}
	  }
	  tred2(V, d, e);
	  tql2(V, d, e);
	}




	DamageTCPlaneStress2DLaw::DamageTCPlaneStress2DLaw()
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
#ifdef DAM_TC_2D_IMPLEX
		, m_rt_converged_old(0.0)
		, m_rc_converged_old(0.0)
		, m_strain()
		, m_strain_converged()
		, m_dTime_n(0.0)
		, m_dTime_n_converged(0.0)
		, m_rt_impl_temp(0.0)
		, m_rc_impl_temp(0.0)
#endif // DAM_TC_2D_IMPLEX
		, m_suggested_time_step(0.0)
		, m_error_code(0.0)
#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
		, m_E(0.0)
		, m_has_changed_reg(false)
		, m_change_reg_t_x(0.0)
		, m_change_reg_c_x(0.0)
#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
	{
	}

	ConstitutiveLaw::Pointer DamageTCPlaneStress2DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer( new DamageTCPlaneStress2DLaw() );
	}

	DamageTCPlaneStress2DLaw::SizeType DamageTCPlaneStress2DLaw::WorkingSpaceDimension()
	{
		return 2;
	}

	DamageTCPlaneStress2DLaw::SizeType DamageTCPlaneStress2DLaw::GetStrainSize()
	{
		return 3;
	}

	bool DamageTCPlaneStress2DLaw::Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == TEMPERATURE)
			return true;
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == DAMAGE_C)
			return true;
		if(rThisVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			return true;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			return true;
		if(rThisVariable == SUGGESTED_TIME_STEP)
			return true;
		return false;
	}

	bool DamageTCPlaneStress2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_STRAIN)
			return true;
		if(rThisVariable == DISCONTINUITY_DIRECTION)
			return true;
		return false;
	}

	bool DamageTCPlaneStress2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool DamageTCPlaneStress2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool DamageTCPlaneStress2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	double& DamageTCPlaneStress2DLaw::GetValue(
		const Variable<double>& rThisVariable,
		double& rValue)
	{
		rValue = 0.0;
		if(rThisVariable == DAMAGE_T || rThisVariable == TEMPERATURE)
			rValue = m_damage_t;
		else if(rThisVariable == DAMAGE_C)
			rValue = m_damage_c;
		else if(rThisVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			rValue = m_lch_multiplier;
		else if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		else if(rThisVariable == SUGGESTED_TIME_STEP)
			rValue = m_suggested_time_step;
		else if(rThisVariable == ENH_STRAIN_PARAM_1)
			rValue = m_rt_converged;
		else if(rThisVariable == ENH_STRAIN_PARAM_2)
			rValue = m_rc_converged;
		return rValue;
	}

	Vector& DamageTCPlaneStress2DLaw::GetValue(
		const Variable<Vector>& rThisVariable,
		Vector& rValue)
	{
		if(rThisVariable == INITIAL_STRAIN) {
			if(rValue.size() != m_initial_strain.size())
				rValue.resize(m_initial_strain.size());
			noalias(rValue) = m_initial_strain;
		}
		else if(rThisVariable == DISCONTINUITY_DIRECTION) {
			rValue = ZeroVector(3);
			if(m_localized) {
				double disc_len = 1.0;
				rValue(0) = disc_len*std::cos(m_localization_angle);
				rValue(1) = disc_len*std::sin(m_localization_angle);
			}
		}
		return rValue;
	}

	Matrix& DamageTCPlaneStress2DLaw::GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & DamageTCPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & DamageTCPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void DamageTCPlaneStress2DLaw::SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == DAMAGE_T)
			m_damage_t = rValue;
		else if(rVariable == DAMAGE_C)
			m_damage_c = rValue;
		else if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER) {
			m_lch_multiplier = rValue;
#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
			if(!m_has_changed_reg) {
				if(m_lch_multiplier != 1.0) {
					m_has_changed_reg = true;
					// find the equivalent strain at which the regularization has changed
					m_change_reg_t_x = m_rt_converged / m_E;
					m_change_reg_c_x = m_rc_converged / m_E;
				}
			}
#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
		}
		else if(rVariable == ENH_STRAIN_PARAM_1)
			m_rt_converged = rValue;
		else if(rVariable == ENH_STRAIN_PARAM_2)
			m_rc_converged = rValue;
	}

	void DamageTCPlaneStress2DLaw::SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == INITIAL_STRAIN) {
			if(rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
	}

	void DamageTCPlaneStress2DLaw::SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTCPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTCPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool DamageTCPlaneStress2DLaw::ValidateInput(const Properties& rMaterialProperties)
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

	DamageTCPlaneStress2DLaw::StrainMeasure DamageTCPlaneStress2DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	DamageTCPlaneStress2DLaw::StressMeasure DamageTCPlaneStress2DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool DamageTCPlaneStress2DLaw::IsIncremental()
	{
		return false;
	}

	void DamageTCPlaneStress2DLaw::InitializeMaterial(
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

#ifdef DAM_TC_2D_IMPLEX
			m_rt_converged_old = m_rt;
			m_rc_converged_old = m_rc;
			m_strain = ZeroVector(this->GetStrainSize());
			m_strain_converged = ZeroVector(this->GetStrainSize());
			m_dTime_n = 0.0;
			m_dTime_n_converged = 0.0;
#endif // DAM_TC_2D_IMPLEX

			m_localized = false;
			m_localized_converged = false;
			m_localization_angle = 0.0;

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
			m_E = rMaterialProperties[YOUNG_MODULUS];
			m_has_changed_reg = false;
			m_change_reg_t_x = 0.0;
			m_change_reg_c_x = 0.0;
#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

		}
	}

	void DamageTCPlaneStress2DLaw::InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTCPlaneStress2DLaw::FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
#ifdef DAM_TC_2D_IMPLEX

		m_rt = m_rt_impl_temp;
		m_rc = m_rc_impl_temp;

		// move from n to n-1
		m_rt_converged_old  = m_rt_converged;
		m_rc_converged_old  = m_rc_converged;
		m_dTime_n_converged = m_dTime_n;
		m_strain_converged = m_strain;

#endif // DAM_TC_2D_IMPLEX

		// save converged values
		m_rt_converged = m_rt;
		m_rc_converged = m_rc;

		m_localized_converged = m_localized;
	}

	void DamageTCPlaneStress2DLaw::InitializeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTCPlaneStress2DLaw::FinalizeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		const ProcessInfo&  pinfo = rValues.GetProcessInfo();
		const GeometryType& geom  = rValues.GetElementGeometry();
		const Properties&   props = rValues.GetMaterialProperties();

		const Vector& strain_vector       = rValues.GetStrainVector();
		Vector&       stress_vector       = rValues.GetStressVector();

		CalculationData data;
		this->InitializeCalculationData(props, geom, rValues.GetShapeFunctionsValues(), pinfo, data);

		this->CalculateMaterialResponseInternal(strain_vector, stress_vector, data);

		if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
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
		}

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
			double h = 1.0E-8;

			// perturbed vectors
			Vector strain_bar(n);
			Vector S1(n);
			Vector S2(n);

			// apply perturbation to each strain component...
			for(size_t j = 0; j < n; j++)
			{
				noalias(strain_bar) = strain_vector;

				strain_bar(j) = strain_vector(j) + h;
				this->CalculateMaterialResponseInternal(strain_bar, S2, data);

				for(size_t i = 0; i < n; i++)
					constitutive_matrix(i,j) = (S2(i) - stress_vector(i))/h;
			}

			// restore internal variables
			m_rt	   = save_rt;
			m_rc	   = save_rc;
			m_damage_t = save_dt;
			m_damage_c = save_dc;

			// localization direction
			m_localized = m_localized_converged;
			double angle = 0.0;
			double det_Qn = min_det_acoustic_tensor_2d(constitutive_matrix,angle);
			if(det_Qn < 0.0) {
				m_localized = true;
				m_localization_angle = angle;
			}
		}
	}

	void DamageTCPlaneStress2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageTCPlaneStress2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{

	}

	void DamageTCPlaneStress2DLaw::ResetMaterial(
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

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2
		m_E = 0.0;
		m_has_changed_reg = false;
		m_change_reg_t_x = 0.0;
		m_change_reg_c_x = 0.0;
#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

	}

	void DamageTCPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set( PLANE_STRESS_LAW );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int DamageTCPlaneStress2DLaw::Check(
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

		if( !rMaterialProperties.Has(VISCOSITY) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: VISCOSITY", "");

		return 0;

		KRATOS_CATCH("");
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponse(const Vector& StrainVector,
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
		Matrix F(IdentityMatrix(2,2));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);
		parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
		this->CalculateMaterialResponseCauchy(parameters);
	}

	void DamageTCPlaneStress2DLaw::InitializeCalculationData(const Properties& props,
															 const GeometryType& geom,
															 const Vector& N,
															 const ProcessInfo& pinfo,
															 CalculationData& data)
	{
		// random imperfections
		double impf = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(geom, N);

		// elasticity
		data.E   = props[YOUNG_MODULUS];//*impf;
		data.nu  = props[POISSON_RATIO];
		this->CalculateElasticityMatrix(data);

		// tension
		data.ft = props[DAMAGE_STRESS_T_0]*impf;
		data.Gt = props[FRACTURE_ENERGY_T]*impf*impf;
		data.m2  = props.Has(DAMAGE_TENSILE_SURFACE_S1) ? props[DAMAGE_TENSILE_SURFACE_S1] : 1.0;
		data.m2  = std::min(std::max(data.m2,0.0),1.0);

		// compression
		data.fc0 = props[DAMAGE_STRESS_C_0]*impf;
		data.fcp = props[DAMAGE_STRESS_C_P]*impf;
		data.fcr = props[DAMAGE_STRESS_C_R]*impf;
		data.ep  = props[DAMAGE_STRAIN_C_P]*impf;
		data.c1  = props.Has(DAMAGE_COMPRESSIVE_LAW_C1) ? props[DAMAGE_COMPRESSIVE_LAW_C1] : 0.65;
		data.c2  = props.Has(DAMAGE_COMPRESSIVE_LAW_C2) ? props[DAMAGE_COMPRESSIVE_LAW_C2] : 0.50;
		data.c3  = props.Has(DAMAGE_COMPRESSIVE_LAW_C3) ? props[DAMAGE_COMPRESSIVE_LAW_C3] : 1.50;
		data.Gc  = props[FRACTURE_ENERGY_C]*impf*impf;
		data.bm  = props[BIAXIAL_COMPRESSION_MULTIPLIER];
		data.m1  = props.Has(SHEAR_COMPRESSION_REDUCTION) ? props[SHEAR_COMPRESSION_REDUCTION] : 0.5;
		data.m1  = std::min(std::max(data.m1,0.0),1.0);

		// effective stress data
		data.S.resize(3,false);
		data.Si.resize(2,false);
		data.ST.resize(3,false);
		data.SC.resize(3,false);
		data.PT.resize(3,3,false);
		data.PC.resize(3,3,false);

		// misc
		data.eta = props[VISCOSITY]*impf;
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
		2 = generalized rankine
		*/
		data.tensile_damage_model = props.Has(DAMAGE_TENSILE_MODEL) ? props[DAMAGE_TENSILE_MODEL] : 0;

		// data for generalized rankine model
		data.grank_a = props.Has(GENRANKINE_SURFACE_PARAM_A) ? props[GENRANKINE_SURFACE_PARAM_A] : 20.0;
		data.grank_b = props.Has(GENRANKINE_SURFACE_PARAM_B) ? props[GENRANKINE_SURFACE_PARAM_B] : 10.0;
		data.grank_c = props.Has(GENRANKINE_SURFACE_PARAM_C) ? props[GENRANKINE_SURFACE_PARAM_C] :  3.0;
	}

	void DamageTCPlaneStress2DLaw::CalculateElasticityMatrix(CalculationData& data)
	{
		if(data.C0.size1() != 3 || data.C0.size2() != 3)
			data.C0.resize(3,3,false);

		double c1 = data.E / (1.0 - data.nu*data.nu);
		double c2 = c1 * data.nu;
		double c3 = c1 * (1.0 - data.nu) / 2.0;

		data.C0(0,0) = c1;		data.C0(0,1) = c2;		data.C0(0,2) = 0.0;
		data.C0(1,0) = c2;		data.C0(1,1) = c1;		data.C0(1,2) = 0.0;
		data.C0(2,0) = 0.0;		data.C0(2,1) = 0.0;		data.C0(2,2) = c3;
	}

	void DamageTCPlaneStress2DLaw::TensionCompressionSplit(CalculationData& data)
	{
		const Vector& S  = data.S;
		Vector&       Si = data.Si;
		Vector&       ST = data.ST;
		Vector&       SC = data.SC;

		Matrix&       PT = data.PT;
		Matrix&       PC = data.PC;
		array_1d<double,2> N1;
		array_1d<double,2> N2;

		if(Si.size() != 2) Si.resize(2,false);
		if(ST.size() != 3) ST.resize(3,false);
		if(SC.size() != 3) SC.resize(3,false);

		ST.clear();
		SC.clear();

		if(std::abs(S(2)) < DAM_TC_PREC)
		{
			if(S(0) > S(1))
			{
				Si(0) = S(0);
				Si(1) = S(1);
				if(Si(0) < 0.0)
					SC(0) = Si(0);
				else
					ST(0) = Si(0);
				if(Si(1) < 0.0)
					SC(1) = Si(1);
				else
					ST(1) = Si(1);

				N1(0)=1.0; N2(0)=0.0;
				N1(1)=0.0; N2(1)=1.0;
			}
			else
			{
				Si(0) = S(1);
				Si(1) = S(0);
				if(Si(0) < 0.0)
					SC(1) = Si(0);
				else
					ST(1) = Si(0);
				if(Si(1) < 0.0)
					SC(0) = Si(1);
				else
					ST(0) = Si(1);

				N1(0)=0.0; N2(0)=1.0;
				N1(1)=1.0; N2(1)=0.0;
			}
		}
		else
		{
			Matrix V(2,2);
			Matrix A(2,2);
			A(0,0) = S(0);
			A(1,1) = S(1);
			A(0,1) = S(2);
			A(1,0) = S(2);
			eigen_decomposition_sym_2x2(A,V,Si);
			double N1x = V(0,0);
			double N1y = V(1,0);
			double N2x = V(0,1);
			double N2y = V(1,1);
			double norm_N1 = std::sqrt(N1x*N1x + N1y*N1y);
			double norm_N2 = std::sqrt(N2x*N2x + N2y*N2y);
			if(norm_N1 > 0.0)
			{
				N1x /= norm_N1;
				N1y /= norm_N1;
			}
			if(norm_N2 > 0.0)
			{
				N2x /= norm_N2;
				N2y /= norm_N2;
			}
			if(Si(0) < Si(1))
			{
				double aux;
				aux = Si(0); Si(0) = Si(1); Si(1) = aux;
				aux = N1x; N1x = N2x; N2x = aux;
				aux = N1y; N1y = N2y; N2y = aux;
			}


			if(Si(0) < 0.0)
			{
				SC(0) = Si(0)*N1x*N1x;
				SC(1) = Si(0)*N1y*N1y;
				SC(2) = Si(0)*N1x*N1y;
			}
			else
			{
				ST(0) = Si(0)*N1x*N1x;
				ST(1) = Si(0)*N1y*N1y;
				ST(2) = Si(0)*N1x*N1y;
			}
			if(Si(1) < 0.0)
			{
				SC(0) += Si(1)*N2x*N2x;
				SC(1) += Si(1)*N2y*N2y;
				SC(2) += Si(1)*N2x*N2y;
			}
			else
			{
				ST(0) += Si(1)*N2x*N2x;
				ST(1) += Si(1)*N2y*N2y;
				ST(2) += Si(1)*N2x*N2y;
			}

			N1(0)=N1x; N2(0)=N2x;
			N1(1)=N1y; N2(1)=N2y;

		}


		array_1d<double,3> p11;
		p11(0) = N1(0)*N1(0);
		p11(1) = N1(1)*N1(1);
		p11(2) = N1(0)*N1(1);
		array_1d<double,3> p22;
		p22(0) = N2(0)*N2(0);
		p22(1) = N2(1)*N2(1);
		p22(2) = N2(0)*N2(1);

		array_1d<double,3> p11_voigt( p11 );
		p11_voigt(2) *= 2.0;
		array_1d<double,3> p22_voigt( p22 );
		p22_voigt(2) *= 2.0;

		PT.clear();
		noalias(PT) += DAM_TC_HVS(Si(0)) * outer_prod(p11,p11_voigt);
		noalias(PT) += DAM_TC_HVS(Si(1)) * outer_prod(p22,p22_voigt);

		noalias(PC) = IdentityMatrix(3,3) - PT;

	}

	void DamageTCPlaneStress2DLaw::TrialEquivalentStressTension(CalculationData& data, double& rt_trial)
	{
		rt_trial = 0.0;
		if(data.Si(0) > 0.0)
		{
			if(data.tensile_damage_model == 0) // LUBLINER
			{
				double fc = data.fcp*data.m2;
				double fb = fc * data.bm;
				double alpha = (fb-fc)/(2.0*fb-fc);
				double ft = data.ft;
				double I1 = data.Si(0) + data.Si(1);
				double S0_S1 = data.Si(0)-data.Si(1);
				double J2 = 1.0/6.0*( S0_S1*S0_S1 + data.Si(0)*data.Si(0) + data.Si(1)*data.Si(1) );
				double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
				double smax = std::max(std::max(data.Si(0), data.Si(1)),0.0);
				rt_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + beta*smax) /fc*ft;
			}
			else if(data.tensile_damage_model == 1) // RANKINE
			{
				rt_trial = std::max(std::max(data.Si(0), data.Si(1)), 0.0);
			}
			else if(data.tensile_damage_model == 2) // GENERALIZED RANKINE
			{
				double fc =  data.fcp;
				double ft =  data.ft;
				double x0 = -data.grank_a*fc;
				double x1 = -data.grank_b*fc;
				double x2 =  0.0;
				double y0 =  0.0;
				double y1 =  data.grank_c*ft;
				double y2 =  ft;
				double sig_max = std::max(data.Si(0),data.Si(1));
				double sig_min = std::min(data.Si(0),data.Si(1));
				if(sig_min >= 0.0 && sig_max >= 0.0) {
					rt_trial = sig_max;
				}
				else {
					rt_trial=bz_2_f(sig_min,sig_max,x0,x1,x2,y0,y1,y2,ft);
				}
			}
			else // NO DAMAGE IN TENSION
			{
				rt_trial = 0.0;
			}
		}
	}

	void DamageTCPlaneStress2DLaw::TrialEquivalentStressCompression(CalculationData& data, double& rc_trial)
	{
		rc_trial = 0.0;
		if(data.Si(1) < 0.0)
		{
			double fc = data.fc0;
			double fb = fc * data.bm;
			double alpha = (fb-fc)/(2.0*fb-fc);
			double ft = data.ft*data.fc0/data.fcp;

			double I1 = data.Si(0) + data.Si(1);
			double S0_S1 = data.Si(0)-data.Si(1);
			double J2 = 1.0/6.0*( S0_S1*S0_S1 + data.Si(0)*data.Si(0) + data.Si(1)*data.Si(1) );

			double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
			double smax = std::max(std::max(data.Si(0), data.Si(1)),0.0);
			rc_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + data.m1*beta*smax);
		}
	}

	void DamageTCPlaneStress2DLaw::CalculateDamageTension(CalculationData& data, double rt, double& dt)
	{
		if(rt <= data.ft)
		{
			dt = 0.0;
		}
		else
		{

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

			REG2_DIS_VEC_T x;
			REG2_DIS_VEC_T y;

			double lch = data.lch;
			double E   = data.E;
			double ft  = data.ft;
			double G   = data.Gt;

			if(m_has_changed_reg)
			{
				// unsing initial lch
				reg2_discrete_t_law(E,G/m_lch,ft,x,y);
				// find the change point
				unsigned int pn = reg2_find_stop_id(x,y,m_change_reg_t_x);

				// insert point
				double insert_y = reg2_interp_dis_curve(x,y,m_change_reg_t_x);
				std::vector<double>::iterator x_it = x.begin();
				std::vector<double>::iterator y_it = y.begin();
				std::advance(x_it,pn);
				std::advance(y_it,pn);
				x.insert(x_it, m_change_reg_t_x);
				y.insert(y_it, insert_y);
				// calculate post stop inelastic strain
				pn++;

				reg2_regularize_from_point(x,y,pn,m_lch,lch);
				// make sure no snap-back takes place! maybe a warning here...
				reg2_correct_snap_back(x);
				// evaluate the hardening value
				double ix = rt/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dt = 1.0-iy/rt;
			}
			else
			{
				reg2_discrete_t_law(E,G/m_lch,ft,x,y);
				double ix = rt/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dt = 1.0-iy/rt;
			}

#else

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

#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

		}
	}

	void DamageTCPlaneStress2DLaw::CalculateDamageCompression(CalculationData& data, double rc, double& dc)
	{
		if(rc <= data.fc0)
		{
			dc = 0.0;
		}
		else
		{

#ifdef DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

			REG2_DIS_VEC_C x;
			REG2_DIS_VEC_C y;

			double lch = data.lch;
			double E   = data.E;
			double G   = data.Gc;

			if(m_has_changed_reg)
			{
				// unsing initial lch
				reg2_discrete_c_law(data, G/m_lch, x, y);
				// find the change point
				unsigned int pn = reg2_find_stop_id(x,y,m_change_reg_c_x);

				// insert point
				double insert_y = reg2_interp_dis_curve(x,y,m_change_reg_c_x);
				std::vector<double>::iterator x_it = x.begin();
				std::vector<double>::iterator y_it = y.begin();
				std::advance(x_it,pn);
				std::advance(y_it,pn);
				x.insert(x_it, m_change_reg_c_x);
				y.insert(y_it, insert_y);
				// calculate post stop inelastic strain
				pn++;

				reg2_regularize_from_point(x,y,pn,m_lch,lch);
				// make sure no snap-back takes place! maybe a warning here...
				reg2_correct_snap_back(x);
				// evaluate the hardening value
				double ix = rc/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dc = 1.0-iy/rc;
			}
			else
			{
				reg2_discrete_c_law(data, G/m_lch, x, y);
				double ix = rc/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dc = 1.0-iy/rc;
			}

#else

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

#endif // DAM_TC_2D_INCREMENTAL_REGULARIZATION_V2

		}
	}

	void DamageTCPlaneStress2DLaw::CalculateMaterialResponseInternal(const Vector& strain_vector,
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

		if(std::abs(data.S(0)) < DAM_TC_PREC) data.S(0) = 0.0;
		if(std::abs(data.S(1)) < DAM_TC_PREC) data.S(1) = 0.0;
		if(std::abs(data.S(2)) < DAM_TC_PREC) data.S(2) = 0.0;

		this->TensionCompressionSplit(data);

		// compute the equivalent stress measures

		double rt_trial;
		double rc_trial;
		this->TrialEquivalentStressTension(data, rt_trial);
		this->TrialEquivalentStressCompression(data, rc_trial);

		// damage update

#ifdef DAM_TC_2D_IMPLEX

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

#endif // DAM_TC_2D_IMPLEX

		// calculation of stress tensor
		noalias(stress_vector)  = (1.0 - m_damage_t)*data.ST;
		noalias(stress_vector) += (1.0 - m_damage_c)*data.SC;
	}


} /* namespace Kratos.*/
