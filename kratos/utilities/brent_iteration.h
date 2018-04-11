//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_BRENT_ITERATION_H_INCLUDED)
#define  KRATOS_BRENT_ITERATION_H_INCLUDED

#include <cmath>
#include <algorithm>
#include "includes/define.h"

namespace Kratos
{

/**
 Find a root of a one-variable function using Brent's algorithm.
 This is a C++ implementation of the algorithm in Press, Teukolsky,
 Vetterling, Flannery. Numerical Recipes in Fortran (2nd Edition).
 */
class BrentIteration
{

public:

    /**
     Find a root of Funtion in the range defined by Guess1 and Guess2.
     Function is something that takes a double and returns a double.
     Note that Guess1 and Guess2 must be on opposite sides of the
     root, that is Function(Guess1)*Function(Guess2) < 0.0.
     */
	template< typename func_t >
	static double FindRoot(
			func_t Function,
			double Guess1,
			double Guess2,
			double Tolerance,
			int MaximumIterations)
	{
		double a = Guess1;
		double b = Guess2;
		double fa = Function(a);
		double fb = Function(b);

		if (fa*fb > 0.0)
		{
			KRATOS_THROW_ERROR(std::runtime_error,"Error in BrentIteration::FindRoot: The images of both initial guesses have the same sign.","");
		}

		double c = b;
		double fc = fb;
		double p,q,r,s,tol1,xm;
		// Explicit initialization to silence a gcc compilation warning
		double d = 0.0, e = 0.0;
		double eps = (Tolerance < 1e-8) ? Tolerance : 1e-8; // Small value, was machine precision in the original algorithm

		for( int iter = 0; iter < MaximumIterations; iter++ )
		{
			// Reorder values and update interval length
			if ( fb*fc > 0.0 )
			{
				c = a;
				fc = fa;
				d = b-a;
				e = d;
			}

			if ( std::fabs(fc) < std::fabs(fb) )
			{
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}

			// Convergence check
			tol1 = 2.*eps*std::fabs(b) + 0.5*Tolerance;
			xm = 0.5*(c-b);

			if (std::fabs(xm) < tol1 || std::fabs(b) < eps)
				return b;

			if (std::fabs(e) > tol1 && std::fabs(fa) > std::fabs(fb) )
			{
				// Try to use inverse quadratic interpolation
				s = fb/fa;
				if ( std::fabs(a-c) < eps )
				{
					p = 2.0*xm*s;
					q = 1.0-s;
				}
				else
				{
					q = fa/fc;
					r = fb/fc;
					p = s*( 2.0*xm*q*(q-r) - (b-a)*(r-1.0) );
					q = (q-1.0)*(r-1.0)*(s-1.0);
				}

				// Change sign of update if necessary
				if ( p > 0.0 ) 
					q = -q;
				p = std::fabs(p);

				double bound1 = 3.0*xm*q - std::fabs(tol1*q);
				double bound2 = std::fabs(e*q);
				if ( 2.0*p < std::min( bound1 , bound2 ) )
				{
					// Accept inverse quadratic interpolation guess
					e = d;
					d = p/q;
				}
				else
				{
					// Fall back to bisection
					d = xm;
					e = d;
				}
			}
			else
			{
				// Bounds decreasing too slowly, switch to bisection
				d = xm;
				e = d;
			}

			// Move last best guess to a
			a = b;
			fa = fb;

			// Evaluate new guess
			if (std::abs(d) > tol1)
			{
				b += d;
			}
			else
			{
				b += (xm > 0.0) ? std::fabs(tol1) : -std::fabs(tol1);
			}
			fb = Function(b);
		}


		return b;
	}

};

}

#endif // KRATOS_BRENT_ITERATION_H_INCLUDED
