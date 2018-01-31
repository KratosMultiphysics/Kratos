//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//


#if !defined(KRATOS_MESHMOVING_UTILITIES_H_INCLUDED )
#define  KRATOS_MESHMOVING_UTILITIES_H_INCLUDED



// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "ale_application.h"

namespace Kratos
{
	namespace MoveMeshUtilities
	{

        typedef Element BaseType;
        typedef BaseType::MatrixType MatrixType;
		typedef Element::GeometryType GeometryType;
        typedef GeometryData::IntegrationMethod IntegrationMethod;

    	//typedef Properties PropertiesType;

        MatrixType CalculateShapeFunctionDerivatives(const int& rdimension, const double& rPointNumber, GeometryType& rGeometry);

        void CalculateInitialJacobian(MatrixType& rJ0, const double& rPointNumber, GeometryType& rGeometry);

	}  // namespace Move Mesh Utilities.
  
}  // namespace Kratos.

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined
