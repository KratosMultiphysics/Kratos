//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi (porting from utility made by Fabio Petracca and Klaus Sautter)
//

#if !defined(KRATOS_SHELL_UTILITIES_H_INCLUDED )
#define  KRATOS_SHELL_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


namespace Kratos
{
	namespace Utilities
	{
		template<class TVec>
		inline void ShapeFunc(double xi, double eta, TVec & N)
		{
			N(0) = 0.25 * (1.0 - xi) * (1.0 - eta); // node 1
			N(1) = 0.25 * (1.0 + xi) * (1.0 - eta); // node 2
			N(2) = 0.25 * (1.0 + xi) * (1.0 + eta); // node 3
			N(3) = 0.25 * (1.0 - xi) * (1.0 + eta); // node 4
		}

		template<class TMat>
		inline void ShapeFunc_NaturalDerivatives(double xi, double eta,
			TMat & dN)
		{
			dN(0, 0) = -(1.0 - eta) * 0.25;
			dN(1, 0) = (1.0 - eta) * 0.25;
			dN(2, 0) = (1.0 + eta) * 0.25;
			dN(3, 0) = -(1.0 + eta) * 0.25;

			dN(0, 1) = -(1.0 - xi)  * 0.25;
			dN(1, 1) = -(1.0 + xi)  * 0.25;
			dN(2, 1) = (1.0 + xi)  * 0.25;
			dN(3, 1) = (1.0 - xi)  * 0.25;
		}

		inline double dN_seren_dxi(const int actualNodeNumber,const double xi,
			const double eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1.0)*(-0.25*xi + 0.25) -
					0.25*(-eta + 1.0)*(-eta - xi - 1.0);
				break;
			case 2:
				returnValue = (-eta + 1.0)*(0.25*xi + 0.25) +
					0.25*(-eta + 1.0)*(-eta + xi - 1.0);
				break;
			case 3:
				returnValue = (eta + 1.0)*(0.25*xi + 0.25) +
					0.25*(eta + 1.0)*(eta + xi - 1.0);
				break;
			case 4:
				returnValue = -(eta + 1.0)*(-0.25*xi + 0.25) -
					0.25*(eta + 1.0)*(eta - xi - 1.0);
				break;
			case 5:
				returnValue = -1.0*xi*(-eta + 1.0);
				break;
			case 6:
				returnValue = -0.5*eta*eta + 0.5;
				break;
			case 7:
				returnValue = -1.0*xi*(eta + 1.0);
				break;
			case 8:
				returnValue = 0.5*eta*eta - 0.5;
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return returnValue;
		}

		inline double dN_seren_deta(const int actualNodeNumber,const double xi,
			const double eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1.0)*(-0.25*xi + 0.25) -
					(-0.25*xi + 0.25)*(-eta - xi - 1.0);
				break;
			case 2:
				returnValue = -(-eta + 1.0)*(0.25*xi + 0.25) -
					(0.25*xi + 0.25)*(-eta + xi - 1.0);
				break;
			case 3:
				returnValue = (eta + 1.0)*(0.25*xi + 0.25) +
					(0.25*xi + 0.25)*(eta + xi - 1.0);
				break;
			case 4:
				returnValue = (eta + 1.0)*(-0.25*xi + 0.25) +
					(-0.25*xi + 0.25)*(eta - xi - 1.0);
				break;
			case 5:
				returnValue = 0.5*xi*xi - 0.5;
				break;
			case 6:
				returnValue = -1.0*eta*(xi + 1.0);
				break;
			case 7:
				returnValue = -0.5*xi*xi + 0.5;
				break;
			case 8:
				returnValue = -1.0*eta*(-xi + 1.0);
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return returnValue;
		}
		
		inline void InterpToStandardGaussPoints(double& v1, double& v2,
			double& v3)
		{
			double vg1 = v1;
			double vg2 = v2;
			double vg3 = v3;
#ifdef OPT_AVERAGE_RESULTS
			v1 = (vg1 + vg2 + vg3) / 3.0;
			v2 = (vg1 + vg2 + vg3) / 3.0;
			v3 = (vg1 + vg2 + vg3) / 3.0;
#else
			v1 = (2.0*vg1) / 3.0 - vg2 / 3.0 + (2.0*vg3) / 3.0;
			v2 = (2.0*vg1) / 3.0 + (2.0*vg2) / 3.0 - vg3 / 3.0;
			v3 = (2.0*vg2) / 3.0 - vg1 / 3.0 + (2.0*vg3) / 3.0;
#endif // OPT_AVERAGE_RESULTS
		}

		inline void InterpToStandardGaussPoints(std::vector< double >& v)
		{
			if (v.size() != 3) return;
			InterpToStandardGaussPoints(v[0], v[1], v[2]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double,
			3> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 3; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double,
			6> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 6; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Vector >& v)
		{
			if (v.size() != 3) return;
			size_t ncomp = v[0].size();
			for (int i = 1; i < 3; i++)
				if (v[i].size() != ncomp)
					return;
			for (size_t i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Matrix >& v)
		{
			if (v.size() != 3) return;
			size_t nrows = v[0].size1();
			size_t ncols = v[0].size2();
			for (int i = 1; i < 3; i++)
				if (v[i].size1() != nrows || v[i].size2() != ncols)
					return;
			for (size_t i = 0; i < nrows; i++)
				for (size_t j = 0; j < ncols; j++)
					InterpToStandardGaussPoints
					(v[0](i, j), v[1](i, j), v[2](i, j));
		}

	}
  
}  // namespace Kratos.

#endif // KRATOS_SHELL_UTILITIES_H_INCLUDED  defined 


