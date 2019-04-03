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
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"
#include "includes/variables.h"

#include "custom_utilities/imperfection_utilities.h"

#include <math.h>

#ifdef MULTISCALE_APPLICATION_USE_EIGEN
#include <external_includes/eigen-3.2.2/Eigen/Dense>
#endif // MULTISCALE_APPLICATION_USE_EIGEN



#define DAM_TC_PREC 0.0

#define DAM_TC_SIGN(X) (X == 0.0 ? 0.0 : ( X > 0.0 ? 1.0 : -1.0 ))
#define DAM_TC_HVS(X)  ( X > 0.0 ? 1.0 : ( X < 0.0 ? -1.0 : 0.5) )
#define DAM_TC_POS(X)  ( X > 0.0 ? X : 0.0 )

#define DAM_TC_OPTIMIZE_LCH

namespace Kratos
{

	// =====================================================================================================
	//
	// UTILITIES
	//
	// =====================================================================================================

	// todo: remove this. only used for very deformed rectangular geometries

#ifdef DAM_TC_OPTIMIZE_LCH

	inline double get_optimized_lch(const Element::GeometryType& geom)
	{
		double lch = geom.Length();
		if(geom.WorkingSpaceDimension() == 3) {
			//if(geom.PointsNumber() == 4 /*&& geom.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral*/) {
			//	// 3D 4N (shell)
			//	double ax = (geom[0].X0()+geom[3].X0())/2.0;
			//	double ay = (geom[0].Y0()+geom[3].Y0())/2.0;
			//	double az = (geom[0].Z0()+geom[3].Z0())/2.0;
			//	double bx = (geom[1].X0()+geom[2].X0())/2.0;
			//	double by = (geom[1].Y0()+geom[2].Y0())/2.0;
			//	double bz = (geom[1].Z0()+geom[2].Z0())/2.0;
			//	double cx = (geom[0].X0()+geom[1].X0())/2.0;
			//	double cy = (geom[0].Y0()+geom[1].Y0())/2.0;
			//	double cz = (geom[0].Z0()+geom[1].Z0())/2.0;
			//	double dx = (geom[2].X0()+geom[3].X0())/2.0;
			//	double dy = (geom[2].Y0()+geom[3].Y0())/2.0;
			//	double dz = (geom[2].Z0()+geom[3].Z0())/2.0;
			//	double v1x = bx-ax;
			//	double v1y = by-ay;
			//	double v1z = bz-az;
			//	double v2x = dx-cx;
			//	double v2y = dy-cy;
			//	double v2z = dz-cz;
			//	double lx = std::sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
			//	double ly = std::sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
			//	lch = std::min(lx,ly);
			//}
			//
			double min_x = std::numeric_limits<double>::max();
			double max_x = -min_x;
			double min_y = min_x;
			double max_y = max_x;
			for(size_t i = 0; i < geom.size(); i++) {
				double ix = geom[i].X0();
				double iy = geom[i].Y0();
				min_x = std::min(min_x, ix);
				max_x = std::max(max_x, ix);
				min_y = std::min(min_y, iy);
				max_y = std::max(max_y, iy);
			}
			double lx = max_x-min_x;
			double ly = max_y-min_y;
			double lch = (lx+ly)/2.0;
		}
		return lch;
	}

#endif // DAM_TC_OPTIMIZE_LCH

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
	inline void dam_tc_pii(const array_1d<double,3>& ni, Matrix& pii)
	{
		double A1 = ni(2)*ni(2);
		double A2 = ni(1)*ni(1);
		double A3 = ni(0)*ni(0);
		double A4 = ni(2);
		double A5 = ni(0);
		double A6 = ni(1);
		double A7 = 2.0*A1*A5*A6;
		double A8 = 2.0*A3*A4*A6;
		double A9 = 2.0*A2*A4*A5;
		double A10 = A1*A2;
		pii(0, 0) = A3*A3;
		pii(0, 1) = A2*A3;
		pii(0, 2) = A1*A3;
		pii(0, 3) = 2.0*A3*A5*A6;
		pii(0, 4) = A8;
		pii(0, 5) = 2.0*A3*A4*A5;
		pii(1, 0) = A2*A3;
		pii(1, 1) = A2*A2;
		pii(1, 2) = A10;
		pii(1, 3) = 2.0*A2*A5*A6;
		pii(1, 4) = 2.0*A2*A4*A6;
		pii(1, 5) = A9;
		pii(2, 0) = A1*A3;
		pii(2, 1) = A10;
		pii(2, 2) = A1*A1;
		pii(2, 3) = A7;
		pii(2, 4) = 2.0*A1*A4*A6;
		pii(2, 5) = 2.0*A1*A4*A5;
		pii(3, 0) = A3*A5*A6;
		pii(3, 1) = A2*A5*A6;
		pii(3, 2) = A1*A5*A6;
		pii(3, 3) = 2.0*A2*A3;
		pii(3, 4) = A9;
		pii(3, 5) = A8;
		pii(4, 0) = A3*A4*A6;
		pii(4, 1) = A2*A4*A6;
		pii(4, 2) = A1*A4*A6;
		pii(4, 3) = A9;
		pii(4, 4) = 2.0*A1*A2;
		pii(4, 5) = A7;
		pii(5, 0) = A3*A4*A5;
		pii(5, 1) = A2*A4*A5;
		pii(5, 2) = A1*A4*A5;
		pii(5, 3) = A8;
		pii(5, 4) = A7;
		pii(5, 5) = 2.0*A1*A3;
	}
	inline void dam_tc_pij(const array_1d<double,3>& ni, const array_1d<double,3>& nj, Matrix& pij)
	{
		double A1 = ni(0)*nj(2) + ni(2)*nj(0);
		double A2 = ni(1)*nj(2) + ni(2)*nj(1);
		double A3 = ni(0)*nj(1) + ni(1)*nj(0);
		double A4 = nj(2);
		double A5 = ni(2);
		double A6 = nj(1);
		double A7 = ni(1);
		double A8 = nj(0);
		double A9 = ni(0);
		double A10 = A4*A4;
		double A11 = A5*A5;
		double A12 = A6*A6;
		double A13 = A7*A7;
		double A14 = A8*A8;
		double A15 = A9*A9;
		double A16 = (A1*A2)/2.0;
		double A17 = (A1*A3)/2.0;
		double A18 = (A2*A3)/2.0;
		double A19 = A4*A5*A6*A7;
		double A20 = A4*A5*A8*A9;
		pij(0, 0) = A14*A15;
		pij(0, 1) = A6*A7*A8*A9;
		pij(0, 2) = A20;
		pij(0, 3) = A6*A8*A15 + A7*A9*A14;
		pij(0, 4) = A2*A8*A9;
		pij(0, 5) = A4*A8*A15 + A5*A9*A14;
		pij(1, 0) = A6*A7*A8*A9;
		pij(1, 1) = A12*A13;
		pij(1, 2) = A19;
		pij(1, 3) = A6*A8*A13 + A7*A9*A12;
		pij(1, 4) = A4*A6*A13 + A5*A7*A12;
		pij(1, 5) = A1*A6*A7;
		pij(2, 0) = A20;
		pij(2, 1) = A19;
		pij(2, 2) = A10*A11;
		pij(2, 3) = A3*A4*A5;
		pij(2, 4) = A4*A6*A11 + A5*A7*A10;
		pij(2, 5) = A4*A8*A11 + A5*A9*A10;
		pij(3, 0) = (A6*A8*A15)/2.0 + (A7*A9*A14)/2.0;
		pij(3, 1) = (A6*A8*A13)/2.0 + (A7*A9*A12)/2.0;
		pij(3, 2) = (A3*A4*A5)/2.0;
		pij(3, 3) = A3*A3/2.0;
		pij(3, 4) = A18;
		pij(3, 5) = A17;
		pij(4, 0) = (A2*A8*A9)/2.0;
		pij(4, 1) = (A4*A6*A13)/2.0 + (A5*A7*A12)/2.0;
		pij(4, 2) = (A4*A6*A11)/2.0 + (A5*A7*A10)/2.0;
		pij(4, 3) = A18;
		pij(4, 4) = A2*A2/2.0;
		pij(4, 5) = A16;
		pij(5, 0) = (A4*A8*A15)/2.0 + (A5*A9*A14)/2.0;
		pij(5, 1) = (A1*A6*A7)/2.0;
		pij(5, 2) = (A4*A8*A11)/2.0 + (A5*A9*A10)/2.0;
		pij(5, 3) = A17;
		pij(5, 4) = A16;
		pij(5, 5) = A1*A1/2.0;
	}

	// todo: move this into a DAMAGE_INCREMENTAL_FRACTURE_ENERGY_REGULARIZATION_H

#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2

#define REG2_DIS_VEC_T std::vector<double>
#define REG2_DIS_VEC_C std::vector<double>

	template<class T>
	inline unsigned int reg2_find_stop_id(const T& x, const T& y, double xstop)
	{
		/*unsigned int id = -1;
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
		return id;*/

		unsigned int n = x.size();
		unsigned int id = n-1;
		for(unsigned int i=1; i<n; i++) {
			if(y[i]<y[i-1]) {
				id = i-1;
				break;
			}
		}

		for(unsigned int i=id; i<n; i++) {
			if(x[i]>=xstop) {
				id = i;
				break;
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

	inline void reg2_discrete_c_law(const DamageTC3DLaw::CalculationData& data, double Gc,
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
							if(x2-x1 > 0.0)
								iy = (ix-x1)/(x2-x1)*(y2-y1)+y1;
							else
								iy = y1;
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
		if(pn == x.size()-1) return;
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
			S = -0.99;
		for(unsigned int i = pn; i < n; i++)
			x[i] += (x[i]-xp)*S;
	}

#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2

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
#ifdef DAM_TC_3D_SEQLIN
		, m_PT_n(6,6)
		, m_PC_n(6,6)
		, m_PT_n_converged(6,6)
		, m_PC_n_converged(6,6)
		, m_has_n_converged_data(false)
#endif // DAM_TC_3D_SEQLIN
		, m_error_code(0.0)
#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
		, m_E(0.0)
		, m_has_changed_reg(false)
		, m_change_reg_t_x(0.0)
		, m_change_reg_c_x(0.0)
#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
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
		else if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER){
			/*m_lch_multiplier = rValue;*/
#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
			if(!m_has_changed_reg) {
				//if(m_lch_multiplier != 1.0) {
				//	m_has_changed_reg = true;
				//	// find the equivalent strain at which the regularization has changed
				//	m_change_reg_t_x = m_rt_converged / m_E;
				//	m_change_reg_c_x = m_rc_converged / m_E;
				//}
				if(  std::abs(m_lch_multiplier-rValue) > 1.0e-10 /*m_lch_multiplier != rValue*/) {
					m_lch_multiplier = rValue;
					m_has_changed_reg = true;
					// find the equivalent strain at which the regularization has changed
					m_change_reg_t_x = m_rt_converged / m_E;
					m_change_reg_c_x = m_rc_converged / m_E;
				}
			}
#else
			m_lch_multiplier = rValue;
#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
		}
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

			//m_lch = 1.0; // occhio

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

#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
			m_E = rMaterialProperties[YOUNG_MODULUS];
			m_has_changed_reg = false;
			m_change_reg_t_x = 0.0;
			m_change_reg_c_x = 0.0;
#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
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
		// save implicit values
		m_rt = m_rt_impl_temp;
		m_rc = m_rc_impl_temp;
		// move from n to n-1
		m_rt_converged_old  = m_rt_converged;
		m_rc_converged_old  = m_rc_converged;
		m_dTime_n_converged = m_dTime_n;
		m_strain_converged  = m_strain;
#endif // DAM_TC_3D_IMPLEX

#ifdef DAM_TC_3D_SEQLIN
		noalias(m_PT_n_converged) = m_PT_n;
		noalias(m_PC_n_converged) = m_PC_n;
#endif // DAM_TC_3D_SEQLIN

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
		this->InitializeCalculationData(props, geom, rValues.GetShapeFunctionsValues(), pinfo, data);
		this->CalculateMaterialResponseInternal(strain_vector, stress_vector, rValues.GetConstitutiveMatrix(), data);
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

#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
		m_E = 0.0;
		m_has_changed_reg = false;
		m_change_reg_t_x = 0.0;
		m_change_reg_c_x = 0.0;
#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
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

		if( !rMaterialProperties.Has(VISCOSITY) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: VISCOSITY", "");

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
		data.Gt = props[FRACTURE_ENERGY_T]*impf*impf;
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
		data.Gc  = props[FRACTURE_ENERGY_C]*impf*impf;
		data.bm  = props[BIAXIAL_COMPRESSION_MULTIPLIER];
		data.m1  = props.Has(SHEAR_COMPRESSION_REDUCTION) ? props[SHEAR_COMPRESSION_REDUCTION] : 0.5;
		data.m1  = std::min(std::max(data.m1,0.0),1.0);

		// lubliner parameter for compressive states
		data.Kc  = props.Has(LUBLINER_SURFACE_PARAM_KC) ? props[LUBLINER_SURFACE_PARAM_KC] : 2.0/3.0;

		// effective stress data
		data.S.resize(6,false);
		data.Si.resize(3,false);
		data.ST.resize(6,false);
		data.SC.resize(6,false);
		data.PT.resize(6,6,false);
		data.PC.resize(6,6,false);
		data.QT.resize(6,6,false);
		data.QC.resize(6,6,false);

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
		Matrix&       PT = data.PT;
		Matrix&       PC = data.PC;
		Matrix&       QT = data.QT;
		Matrix&       QC = data.QC;

		ST.clear();
		SC.clear();
		PT.clear();
		PC.clear();
		QT.clear();
		QC.clear();

		array_1d<double,3>& v0 = data.v0;
		array_1d<double,3>& v1 = data.v1;
		array_1d<double,3>& v2 = data.v2;
		array_1d<double,3> vtemp;

		// compute eigenvalues and eigenvectors
		double p1 = S(3)*S(3) + S(4)*S(4) + S(5)*S(5);
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
			}

			v0 /= norm_2(v0);
			v1 /= norm_2(v1);
			v2 /= norm_2(v2);
		}

		// order eigenvalues and eigenvectors : 1>2>3
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

		// construct matrices PT and PC
		Matrix pii(6,6);
		// p11
		dam_tc_pii(v0, pii);
		double hs1 = DAM_TC_HVS(Si(0));
		if(hs1 > 0.0) PT += hs1*pii;
		// p22
		dam_tc_pii(v1, pii);
		double hs2 = DAM_TC_HVS(Si(1));
		if(hs2 > 0.0) PT += hs2*pii;
		// p33
		dam_tc_pii(v2, pii);
		double hs3 = DAM_TC_HVS(Si(2));
		if(hs3 > 0.0) PT += hs3*pii;
		// PC
		noalias(PC) = IdentityMatrix(6,6) - PT;

		// construct ST and SC
		noalias(ST) = prod( PT, S );
		noalias(SC) = prod( PC, S );

		// construct QT and QC
		// QT
		noalias(QT) = PT;
		// p12
		dam_tc_pij(v0, v1, pii);
		double s12 = Si(0) - Si(1);
		double pos_s12 = DAM_TC_POS(Si(0)) - DAM_TC_POS(Si(1));
		if(std::abs(s12) > 0.0)
			QT += 2.0*pos_s12/s12*pii;
		else
			QT += pii;
		// p13
		dam_tc_pij(v0, v2, pii);
		double s13 = Si(0) - Si(2);
		double pos_s13 = DAM_TC_POS(Si(0)) - DAM_TC_POS(Si(2));
		if(std::abs(s13) > 0.0)
			QT += 2.0*pos_s13/s13*pii;
		else
			QT += pii;
		// p23
		dam_tc_pij(v1, v2, pii);
		double s23 = Si(1) - Si(2);
		double pos_s23 = DAM_TC_POS(Si(1)) - DAM_TC_POS(Si(2));
		if(std::abs(s23) > 0.0)
			QT += 2.0*pos_s23/s23*pii;
		else
			QT += pii;
		// QC
		noalias(QC) = IdentityMatrix(6,6) - QT;
	}

	void DamageTC3DLaw::TrialEquivalentStressTension(CalculationData& data, double& rt_trial, Vector& d_rt_trial_d_sigma)
	{
		rt_trial = 0.0;
		if(d_rt_trial_d_sigma.size() != 6) d_rt_trial_d_sigma.resize(6, false);
		d_rt_trial_d_sigma.clear();

		if(data.Si(0) > 0.0)
		{
			if(data.tensile_damage_model == 0) // LUBLINER
			{
				// equivalent stress
				double fc    = data.fcp*data.m2;
				double fb    = fc * data.bm;
				double alpha = (fb-fc)/(2.0*fb-fc);
				double ft    = data.ft;
				double gamma = 3.0*(1.0-data.Kc)/(2.0*data.Kc-1.0);
				double I1 = data.Si(0) + data.Si(1) + data.Si(2);
				array_1d<double,3>Sdev;
				for(size_t i=0; i < 3; i++)
					Sdev(i) = data.Si(i)-I1/3.0;
				double J2 = 0.5*inner_prod(Sdev,Sdev);
				double beta = fc/ft*(1.0-alpha)-(1.0+alpha);
				double smax_alg = data.Si(0);
				double smax = DAM_TC_POS( smax_alg);
				double smin = DAM_TC_POS(-smax_alg);
				rt_trial = 1.0/(1.0-alpha)*(alpha*I1 + std::sqrt(3.0*J2) + beta*smax - gamma*smin) /fc*ft;
				// derivative
				// d_rt_trial_d_sigma = alpha*I + beta*H(s0)*(v0 x v0) - gamma*(-H(-s0))*(v0 x v0);
				array_1d<double,3>& v0 = data.v0;
				double Hsmax = DAM_TC_HVS( smax_alg);
				double Hsmin = DAM_TC_HVS(-smax_alg);
				double coeff_v0_tens_v0 = beta*Hsmax + gamma*Hsmin;
				d_rt_trial_d_sigma(0) = alpha +     v0(0)*v0(0)*coeff_v0_tens_v0;
				d_rt_trial_d_sigma(1) = alpha +     v0(1)*v0(1)*coeff_v0_tens_v0;
				d_rt_trial_d_sigma(2) = alpha +     v0(2)*v0(2)*coeff_v0_tens_v0;
				d_rt_trial_d_sigma(3) =         2.0*v0(0)*v0(1)*coeff_v0_tens_v0;
				d_rt_trial_d_sigma(4) =         2.0*v0(1)*v0(2)*coeff_v0_tens_v0;
				d_rt_trial_d_sigma(5) =         2.0*v0(0)*v0(2)*coeff_v0_tens_v0;
				if(std::abs(J2) > 0.0)
				{
					// d_rt_trial_d_sigma += 3/(2*sqrt(3*J2))*SDev.*[1 1 1 2 2 2];
					Vector& S = data.S;
					double d_q_d_j2 = 3.0/(2.0*std::sqrt(3.0*J2));
					d_rt_trial_d_sigma(0) += d_q_d_j2 * (S(0) - I1/3.0);
					d_rt_trial_d_sigma(1) += d_q_d_j2 * (S(1) - I1/3.0);
					d_rt_trial_d_sigma(2) += d_q_d_j2 * (S(2) - I1/3.0);
					d_rt_trial_d_sigma(3) += d_q_d_j2 * (2.0*S(3));
					d_rt_trial_d_sigma(4) += d_q_d_j2 * (2.0*S(4));
					d_rt_trial_d_sigma(5) += d_q_d_j2 * (2.0*S(5));
				}
				d_rt_trial_d_sigma *= 1.0/(1.0-alpha)/fc*ft;
			}
			else if(data.tensile_damage_model == 1) // RANKINE
			{
				// equivalent stress
				rt_trial = data.Si(0);
				// derivative = v0 tens v0
				array_1d<double,3>& v0 = data.v0;
				d_rt_trial_d_sigma(0) = v0(0)*v0(0);
				d_rt_trial_d_sigma(1) = v0(1)*v0(1);
				d_rt_trial_d_sigma(2) = v0(2)*v0(2);
				d_rt_trial_d_sigma(3) = 2.0*v0(0)*v0(1);
				d_rt_trial_d_sigma(4) = 2.0*v0(1)*v0(2);
				d_rt_trial_d_sigma(5) = 2.0*v0(0)*v0(2);
			}
		}
	}

	void DamageTC3DLaw::TrialEquivalentStressCompression(CalculationData& data, double& rc_trial, Vector& d_rc_trial_d_sigma)
	{
		rc_trial = 0.0;
		if(d_rc_trial_d_sigma.size() != 6) d_rc_trial_d_sigma.resize(6, false);
		d_rc_trial_d_sigma.clear();

		if(data.Si(2) < 0.0)
		{
			// equivalent stress
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
			// derivative
			// d_rc_trial_d_sigma = alpha*I + m1*beta*H(s0)*(v0 x v0) - gamma*(-H(-s0))*(v0 x v0);
			array_1d<double,3>& v0 = data.v0;
			double Hsmax = DAM_TC_HVS( smax_alg);
			double Hsmin = DAM_TC_HVS(-smax_alg);
			double coeff_v0_tens_v0 = data.m1*beta*Hsmax + gamma*Hsmin;
			d_rc_trial_d_sigma(0) = alpha +     v0(0)*v0(0)*coeff_v0_tens_v0;
			d_rc_trial_d_sigma(1) = alpha +     v0(1)*v0(1)*coeff_v0_tens_v0;
			d_rc_trial_d_sigma(2) = alpha +     v0(2)*v0(2)*coeff_v0_tens_v0;
			d_rc_trial_d_sigma(3) =         2.0*v0(0)*v0(1)*coeff_v0_tens_v0;
			d_rc_trial_d_sigma(4) =         2.0*v0(1)*v0(2)*coeff_v0_tens_v0;
			d_rc_trial_d_sigma(5) =         2.0*v0(0)*v0(2)*coeff_v0_tens_v0;
			if(std::abs(J2) > 0.0)
			{
				// d_rt_trial_d_sigma += 3/(2*sqrt(3*J2))*SDev.*[1 1 1 2 2 2];
				Vector& S = data.S;
				double d_q_d_j2 = 3.0/(2.0*std::sqrt(3.0*J2));
				d_rc_trial_d_sigma(0) += d_q_d_j2 * (S(0) - I1/3.0);
				d_rc_trial_d_sigma(1) += d_q_d_j2 * (S(1) - I1/3.0);
				d_rc_trial_d_sigma(2) += d_q_d_j2 * (S(2) - I1/3.0);
				d_rc_trial_d_sigma(3) += d_q_d_j2 * (2.0*S(3));
				d_rc_trial_d_sigma(4) += d_q_d_j2 * (2.0*S(4));
				d_rc_trial_d_sigma(5) += d_q_d_j2 * (2.0*S(5));
			}
			d_rc_trial_d_sigma *= 1.0/(1.0-alpha)/fc*ft;
		}
	}

	void DamageTC3DLaw::CalculateDamageTension(CalculationData& data, double rt, double& dt, double& d_dt_d_rt)
	{
		if(rt <= data.ft)
		{
			dt = 0.0;
			d_dt_d_rt = 0.0;
		}
		else
		{
#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2

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
				//pn++;
				//if(pn > x.size()-1) pn = x.size()-1;

				reg2_regularize_from_point(x,y,pn,m_lch,lch);
				// make sure no snap-back takes place! maybe a warning here...
				//reg2_correct_snap_back(x);
				// evaluate the hardening value
				double ix = rt/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dt = 1.0-iy/rt;
				// derivative
				double pert = 1.0e-6*ft;
				double rt_pert = rt+pert;
				double ix_pert = rt_pert/E;
				double iy_pert = reg2_interp_dis_curve(x,y,ix_pert);
				double dt_pert = 1.0-iy_pert/rt_pert;
				d_dt_d_rt = (dt_pert-dt)/pert;
			}
			else
			{
				reg2_discrete_t_law(E,G/m_lch,ft,x,y);
				double ix = rt/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dt = 1.0-iy/rt;
				// derivative
				double pert = 1.0e-6*ft;
				double rt_pert = rt+pert;
				double ix_pert = rt_pert/E;
				double iy_pert = reg2_interp_dis_curve(x,y,ix_pert);
				double dt_pert = 1.0-iy_pert/rt_pert;
				d_dt_d_rt = (dt_pert-dt)/pert;
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
			d_dt_d_rt  = (r0 + A*rt)/(rt*rt*std::exp(A*(rt/r0 - 1.0)));
			double rmin=1.0e-2*ft;
			if((1.0-dt)*rt<rmin) {
				dt = 1.0-rmin/rt;
				d_dt_d_rt = rmin/rt/rt;
			}

			//double e0 = ft/E;
			//double s0 = ft*0.5;
			//double su = 0.0;
			////double eu;
			//double G1 = ft*e0/2.0;
			////double G2 = s0*(eu-e0)/2.0;
			//// G1+G2 = G/lch
			///*
			//G2 = G/lch - G1
			//s0*eu/2.0 -s0*e0/2.0 = G/lch - G1
			//s0*eu/2.0 = G/lch - G1 + s0*e0/2.0
			//*/
			//double eu = (G/lch - G1 + s0*e0/2.0)*2.0/s0;
			//if(rt/E < eu) {
			//	//double ss = s0 - s0/(eu-e0)*(rt/E-e0);
			//	double c1 = s0/(eu-e0)/E;
			//	double c2 = s0/(eu-e0)*e0;
			//	double ss = s0 - rt*c1 + c2;
			//	// ss = (1-dt)*rt;
			//	// s0 -rt*c1 + c2 = rt-d*rt;
			//	// s0 -rt*c1 +c2 -rt = -dt*rt;
			//	// s0 +c2 -rt*(c1+1) = -dt*rt;
			//	// rt*(c1+1) - s0 - c2 = dt*rt;
			//	// dt = (c1+1) -s0/rt - c2/rt
			//	// ss = rt - dt*rt;
			//	// ss -rt = -dt*rt;
			//	// 1 - ss/rt = dt
			//	// dt = 1- ss/rt
			//	dt = 1.0 - ss / rt;
			//	d_dt_d_rt = 0.0;
			//}
			//else {
			//	dt = 1.0;
			//	d_dt_d_rt = 0.0;
			//}
			//double rmin=1.0e-2*ft;
			//if((1.0-dt)*rt<rmin) {
			//	dt = 1.0-rmin/rt;
			//	d_dt_d_rt = rmin/rt/rt;
			//}

#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2

			//G = G/lch;
			//double A1 = ft*ft/2.0/E/G; // 1.0/2.0;
			//double A2 = 0.25; //1.0/6.0;
			//double G1 = A1*G;
			//double G2 = G-G1;
			//double s0 = ft;
			//double e0 = s0/E;
			//double s1 = s0*A2;
			//double e1 = 2.0/s0*(G1 + e0*s1/2.0);
			//if(e1 < e0) e1 = 1.01*e0;
			//double e2 = 2.0*G2/s1;
			//if(e2 < e1) e1 = 1.01*e1;
			//
			//double ix = rt/E;


			//double y      = 0.0;
			//double dy_dix = 0.0;
			//if(ix<=e0) {
			//	y = E*ix;
			//	dy_dix = E;
			//}
			//else if(ix <= e1) {
			//	y = s0+(s1-s0)/(e1-e0)*(ix-e0);
			//	dy_dix = (s0-s1)/(e0-e1);
			//}
			//else if(ix <= e2) {
			//	y = s1-s1/(e2-e1)*(ix-e1);
			//	dy_dix = s1/(e1-e2);
			//}
			//dt        = 1.0-y/rt;
			//d_dt_d_rt = -(dy_dix/E/rt - y/(rt*rt));

			//double rmin=1.0e-2*ft;
			//if((1.0-dt)*rt<rmin) {
			//	dt = 1.0-rmin/rt;
			//	d_dt_d_rt = rmin/rt/rt;
			//}
		}
	}

	void DamageTC3DLaw::CalculateDamageCompression(CalculationData& data, double rc, double& dc, double& d_dc_d_rc)
	{
		if(rc <= data.fc0)
		{
			dc = 0.0;
		}
		else
		{
#ifdef DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2

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
				//pn++;
				if(pn > x.size()-1) pn = x.size()-1;

				reg2_regularize_from_point(x,y,pn,m_lch,lch);
				// make sure no snap-back takes place! maybe a warning here...
				//reg2_correct_snap_back(x);
				// evaluate the hardening value
				double ix = rc/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dc = 1.0-iy/rc;
				// derivative
				double pert = 1.0e-6*data.fcp;
				double rc_pert = rc+pert;
				double ix_pert = rc_pert/E;
				double iy_pert = reg2_interp_dis_curve(x,y,ix_pert);
				double dc_pert = 1.0-iy_pert/rc_pert;
				d_dc_d_rc = (dc_pert-dc)/pert;
			}
			else
			{
				reg2_discrete_c_law(data, G/m_lch, x, y);
				double ix = rc/E;
				double iy = reg2_interp_dis_curve(x,y,ix);
				dc = 1.0-iy/rc;
				// derivative
				double pert = 1.0e-6*data.fcp;
				double rc_pert = rc+pert;
				double ix_pert = rc_pert/E;
				double iy_pert = reg2_interp_dis_curve(x,y,ix_pert);
				double dc_pert = 1.0-iy_pert/rc_pert;
				d_dc_d_rc = (dc_pert-dc)/pert;
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

			// derivative
			double h = (sp-s0)/100.0;
			xi = (rc+h)/E;
			if(xi <= ep)
				s = b3_eval_bezier(xi,e0,sp/E,ep,s0,sp,sp);
			else if(xi <= ek)
				s = b3_eval_bezier(xi,ep,ej,ek,sp,sp,sk);
			else if(xi <= eu)
				s = b3_eval_bezier(xi,ek,er,eu,sk,sr,sr);
			else
				s = sr;
			double dch = 1.0-s/(rc+h);
			d_dc_d_rc = (dch-dc)/h;

#endif // DAM_TC_3D_INCREMENTAL_REGULARIZATION_V2
		}
	}

	void DamageTC3DLaw::CalculateMaterialResponseInternal(const Vector& strain_vector,
														  Vector& stress_vector,
														  Matrix& tangent_matrix,
														  CalculationData& data)
	{
		m_error_code = 0.0;

		size_t strain_size = this->GetStrainSize();

		if(stress_vector.size() != strain_size)
			stress_vector.resize(strain_size,false);

		if(tangent_matrix.size1() != strain_size || tangent_matrix.size2() != strain_size)
			tangent_matrix.resize(strain_size, strain_size, false);

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
		Vector d_rt_trial_d_sigma;
		Vector d_rc_trial_d_sigma;
		this->TrialEquivalentStressTension(data, rt_trial, d_rt_trial_d_sigma);
		this->TrialEquivalentStressCompression(data, rc_trial, d_rc_trial_d_sigma);

		// damage update

		double d_dt_d_rt = 0.0;
		double d_dc_d_rc = 0.0;

#ifdef DAM_TC_3D_IMPLEX

		// time factor
		double time_factor = 0.0;
		if(m_dTime_n_converged>0.0)
			time_factor = data.dTime/m_dTime_n_converged;
		m_dTime_n = data.dTime;

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
		this->CalculateDamageTension(data, m_rt, m_damage_t, d_dt_d_rt);
		this->CalculateDamageCompression(data, m_rc, m_damage_c, d_dc_d_rc);
#else

		if(rt_trial > m_rt)
			m_rt = rate_coeff_1*m_rt + rate_coeff_2*rt_trial;
		this->CalculateDamageTension(data, m_rt, m_damage_t, d_dt_d_rt);

		if(rc_trial > m_rc)
			m_rc = rate_coeff_1*m_rc + rate_coeff_2*rc_trial;
		this->CalculateDamageCompression(data, m_rc, m_damage_c, d_dc_d_rc);

#endif // DAM_TC_3D_IMPLEX

#ifdef DAM_TC_3D_SEQLIN
		noalias(m_PT_n) = data.PT;
		noalias(m_PC_n) = data.PC;
		if(!m_has_n_converged_data)
		{
			noalias(m_PT_n_converged) = m_PT_n;
			noalias(m_PC_n_converged) = m_PC_n;
			m_has_n_converged_data = true;
		}

#ifndef DAM_TC_3D_IMPLEX

		this->CalculateDamageTension(data, m_rt_converged, m_damage_t, d_dt_d_rt);
		this->CalculateDamageCompression(data, m_rc_converged, m_damage_c, d_dc_d_rc);

#endif // !DAM_TC_3D_IMPLEX

		// error control
#ifdef DAM_TC_3D_IMPLEX
		/*double aux_dt = 0.0;
		double aux_dc = 0.0;
		double dummy;
		this->CalculateDamageTension(data, m_rt_converged, aux_dt, dummy);
		this->CalculateDamageCompression(data, m_rc_converged, aux_dc, dummy);
		aux_dt += m_damage_t;
		aux_dc += m_damage_c;
		if(data.dTime >= 1.0e-4)
		{
			double max_dam_incr = 0.1;
			if(aux_dt > max_dam_incr)
				m_error_code = -1.0;
			if(aux_dc > max_dam_incr)
				m_error_code = -1.0;
		}*/
#endif // DAM_TC_3D_IMPLEX

#endif // DAM_TC_3D_SEQLIN

		// calculation of stress tensor
#ifdef DAM_TC_3D_SEQLIN
		noalias(data.ST) = prod(m_PT_n_converged, data.S);
		noalias(data.SC) = prod(m_PC_n_converged, data.S);
#endif // DAM_TC_3D_SEQLIN
		noalias(stress_vector)  = (1.0 - m_damage_t)*data.ST;
		noalias(stress_vector) += (1.0 - m_damage_c)*data.SC;

		// calculation of tangent matrix
		Matrix W( IdentityMatrix(6,6) );
#ifdef DAM_TC_3D_SEQLIN
		if(m_damage_t > 0.0)
			noalias(W) -= m_damage_t*m_PT_n_converged;
		if(m_damage_c > 0.0)
			noalias(W) -= m_damage_c*m_PC_n_converged;
#else
		if(m_damage_t > 0.0)
			noalias(W) -= m_damage_t*data.QT;
		if(m_damage_c > 0.0)
			noalias(W) -= m_damage_c*data.QC;
#ifndef DAM_TC_3D_IMPLEX
		if(d_dt_d_rt != 0.0)
			noalias(W) -= rate_coeff_2*d_dt_d_rt*outer_prod( data.ST, d_rt_trial_d_sigma );
		if(d_dc_d_rc != 0.0)
			noalias(W) -= rate_coeff_2*d_dc_d_rc*outer_prod( data.SC, d_rc_trial_d_sigma );
#endif // !DAM_TC_3D_IMPLEX
#endif // DAM_TC_3D_SEQLIN


		noalias(tangent_matrix) = prod( W, data.C0 );
	}


} /* namespace Kratos.*/
