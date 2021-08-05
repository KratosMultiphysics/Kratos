// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

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
                    const double *xmax, const double *xold1, const double *xold2, const int iter, double *low, double *upp);

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
                    const double *xmax, const double *xold1, const double *xold2, const int iter, double *low, double *upp);

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
// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

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
                    const double *xmax);

        void Reset() { iter = 0; };

        int nano, m, iter;

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