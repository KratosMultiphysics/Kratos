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

#if !defined(KRATOS_MMA_SOLVER_H_INCLUDED)
#define  KRATOS_MMA_SOLVER_H_INCLUDED


#pragma once

#include <vector>
// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

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

class MMASolver
{

    public:
        KRATOS_CLASS_POINTER_DEFINITION(MMASolver);

        /// Default constructor.
        MMASolver(int n, int m, double a = 0.0, double c = 1000.0, double d = 0.0);

        void SetAsymptotes(double init, double decrease, double increase);

        void ConstraintModification(bool conMod) {}

        void Update(double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
                    const double *xmax, double *xold1, double *xold2, const int iter, double *low, double *upp);

        ///void Reset() { iter = 0; };

        int nano, m; ///iter;

        const double xmamieps;
        const double epsimin;

        const double raa0;
        const double move, albefa;
        double asyminit, asymdec, asyminc;

        std::vector<double> a, c, d;
        std::vector<double> y;
        double z;

        std::vector<double> lam, mu, s;
        std::vector<double>  alpha, beta, p0, q0, pij, qij, b, grad, hess;///low, upp,

        ///std::vector<double> xold1, xold2;

        void GenSub(const double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
                    const double *xmax, double *xold1,  double *xold2, const int iter, double *low, double *upp);

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

#endif	/* KRATOS_MMA_SOLVER_H_INCLUDED */