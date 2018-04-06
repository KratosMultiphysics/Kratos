//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "shell_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
    namespace ShellUtilities
	{
        typedef Properties PropertiesType; // TODO remove this?

		double dN_seren_dxi(const int actualNodeNumber,const double xi,
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

		double dN_seren_deta(const int actualNodeNumber,const double xi,
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

		void InterpToStandardGaussPoints(double& v1, double& v2,
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

		void InterpToStandardGaussPoints(std::vector< double >& v)
		{
			if (v.size() != 3) return;
			InterpToStandardGaussPoints(v[0], v[1], v[2]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			3> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 3; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			6> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 6; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Vector >& v)
		{
			if (v.size() != 3) return;
			size_t ncomp = v[0].size();
			for (int i = 1; i < 3; i++)
				if (v[i].size() != ncomp)
					return;
			for (size_t i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Matrix >& v)
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
    } // namespace ShellUtilities
}  // namespace Kratos.


