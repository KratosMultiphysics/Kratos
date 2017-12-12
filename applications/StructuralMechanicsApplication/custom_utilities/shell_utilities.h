//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi (porting from utility made by Fabio Petracca and Peter Wilson)
//					 Philipp Bucher
//

#if !defined(KRATOS_SHELL_UTILITIES_H_INCLUDED )
#define  KRATOS_SHELL_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/properties.h"


namespace Kratos
{
	namespace ShellUtilities
	{
		typedef Element::GeometryType GeometryType;

    	typedef Properties PropertiesType;

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

		double dN_seren_dxi(const int actualNodeNumber, const double xi,
			const double eta);
			
		double dN_seren_deta(const int actualNodeNumber, const double xi,
			const double eta);
		
		void InterpToStandardGaussPoints(double& v1, double& v2,
			double& v3);

		void InterpToStandardGaussPoints(std::vector< double >& v);

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			3> >& v);

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			6> >& v);

		void InterpToStandardGaussPoints(std::vector< Vector >& v);

		void InterpToStandardGaussPoints(std::vector< Matrix >& v);

		void CheckVariables();

    	void CheckDofs(GeometryType& rGeom);

    	void CheckProperties(const Element* pTheElement, const ProcessInfo& rCurrentProcessInfo, 
                             const bool IsThickShell = false);

		void CheckSpecificProperties(const Element* pTheElement, const PropertiesType & rProps, const bool IsThickShell);

	}  // namespace Shell Utilities.
  
}  // namespace Kratos.

#endif // KRATOS_SHELL_UTILITIES_H_INCLUDED  defined 


