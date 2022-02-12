////////////////////////////////////////////////////////////////////////////////
// Copyright © 2018 Jérémie Dumas
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
//
// BETA VERSION  0.99
//
// MMA solver using a dual interior point method
//
// Original code by Niels Aage, February 2013
// Modified to use OpenMP by Jun Wu, April 2017
// Various modifications by Jérémie Dumas, June 2017
//
// The class solves a general non-linear programming problem
// on standard from, i.e. non-linear objective f, m non-linear
// inequality constraints g and box constraints on the n
// design variables xmin, xmax.
//
//        min_x^n f(x)
//        s.t. g_j(x) < 0,   j = 1,m
//        xmin < x_i < xmax, i = 1,n
//
// Each call to Update() sets up and solve the following
// convex subproblem:
//
//   min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)
//
//   s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
//             Lj < alphaj <=  xj <= betaj < Uj,  j = 1,n
//             yi >= 0, i = 1,m
//             z >= 0.
//
// NOTE: a0 == 1 in this implementation !!!!
//
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>

class MMASolver {

  public:
	MMASolver(int n, int m, double a = 0.0, double c = 1000.0, double d = 0.0);

	void SetAsymptotes(double init, double decrease, double increase);

	void ConstraintModification(bool conMod) {}

	void Update(double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
	            const double *xmax);

	void Reset() { iter = 0; };

  private:
	int n, m, iter;

	const double xmamieps;
	const double epsimin;

	const double raa0;
	const double move, albefa;
	double asyminit, asymdec, asyminc;

	std::vector<double> a, c, d;
	std::vector<double> y;
	double z;

	std::vector<double> lam, mu, s;
	std::vector<double> low, upp, alpha, beta, p0, q0, pij, qij, b, grad, hess;

	std::vector<double> xold1, xold2;

	void GenSub(const double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
	            const double *xmax);

	void SolveDSA(double *x);
	void SolveDIP(double *x);

	void XYZofLAMBDA(double *x);

	void DualGrad(double *x);
	void DualHess(double *x);
	void DualLineSearch();
	double DualResidual(double *x, double epsi);

	static void Factorize(double *K, int n);
	static void Solve(double *K, double *x, int n);
};
