// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Philipp Hofer, https://github.com/PhiHo-eng
//                   Erich Wehrle, https://github.com/e-dub
//
// ==============================================================================
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
// GCMMA solver using a dual interior point method
//
// Original MMA code by Niels Aage, February 2013
// Extension to GCMMA by Jérémie Dumas, June 2017
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
#if !defined(KRATOS_GCMMA_SOLVER_H_INCLUDED)
#define  KRATOS_GCMMA_SOLVER_H_INCLUDED


#pragma once

#include <vector>
// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>

// External includes
#include <pybind11/pybind11.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "includes/variables.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{

	///@}
	///@name Type Definitions
	///@{


	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Solution utility to filter results.
	/** Detail class definition.

	*/


	class GCMMASolver 
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION(GCMMASolver);

		GCMMASolver(int n,int m, double a = 0.0, double c = 1000.0, double d = 0.0);

		void SetAsymptotes(double init, double decrease, double increase);

		// Compute [L, U, raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
		double OuterUpdate(double *xmma, const double *xval, double f0x, const double *df0dx,
			const double *fx, const double *dfdx, const double *xmin, const double *xmax, 
			double *xold1, double *xold2, const int iter, double *low, double *upp, double f0app);

		// Update [raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
		double InnerUpdate(double *xmma, double f0xnew, const double *fxnew,
			const double *xval, double f0x, const double *df0dx, const double *fx,
			const double *dfdx, const double *xmin, const double *xmax, double *xold1, double *xold2, 
			const int iter, double *low, double *upp, double f0app);

		// Check whether the new solution is conservative
		bool ConCheck(double f0xnew, const double *fxnew, double f0app) const;

		void Reset() { outeriter = 0; };

	private:
		int nano, m;

		const double raa0eps;
		const double raaeps;
		const double xmamieps;
		const double epsimin;

		const double move, albefa;
		double asyminit, asymdec, asyminc;

		double raa0;
		std::vector<double> raa;

		std::vector<double> a, c, d;
		std::vector<double> y;
		double z;

		std::vector<double> lam, mu, s;
		std::vector<double> alpha, beta, p0, q0, pij, qij, b, grad, hess; 

		double r0;
		std::vector<double> r, fapp;

	private:
		// Compute [low, upp, raa0, raa]
		void Asymp(const double *xval, const double *df0dx,
			const double *dfdx, const double *xmin, const double *xmax, 
			double *xold1, double *xold2, const int iter,  double *low, double *upp);

		// Update [raa0, raa]
		void RaaUpdate(const double *xmma, const double *xval, double f0xnew,
			const double *fxnew, const double *xmin, const double *xmax,  double *low, double *upp, const double f0app);

		// Build CGMMA subproblem
		void GenSub(const double *xval, double f0x, const double *df0dx, const double *fx,
			const double *dfdx, const double *xmin, const double *xmax, 
			double *xold1, double *xold2, const int iter,  double *low, double *upp);

		// Compute [f0app, fapp]
		double ComputeApprox(const double *xmma, double *low, double *upp, double f0app);

		void SolveDSA(double *x, double *low, double *upp);
        void SolveDIP(double *x, double *low, double *upp);

        void XYZofLAMBDA(double *x, double *low, double *upp);

        void DualGrad(double *x, double *low, double *upp);
        void DualHess(double *x, double *low, double *upp);
		void DualLineSearch();
		double DualResidual(double *x, double epsi, double *low, double *upp);

		static void Factorize(double *K, int n);
		static void Solve(double *K, double *x, int n);
 		///@name Operators
        ///@{


        ///@}
        ///@name Operations


            ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{



        ///@{


        ///@}

        protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}
        private:
        ///@name Static Member Variables
        ///@{

        

        
        ///@}
        ///@name Member Variables
        ///@{
            


        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        //FilterFunction& operator=(FilterFunction const& rOther);

        /// Copy constructor.
        //FilterFunction(FilterFunction const& rOther);


        ///@}

}; // Class TopologyFilteringUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_GCMMA_SOLVER_H_INCLUDED */
