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


// System includes


// External includes 


// Project includes
#include "move_mesh_utilities.h"


namespace Kratos
{
    namespace MoveMeshUtilities
	{

    
    void CalculateInitialJacobian(MatrixType &rJ0, const double &rPointNumber, GeometryType& rGeometry)
    {
        KRATOS_TRY;

        const IntegrationMethod ThisIntegrationMethod = rGeometry.GetDefaultIntegrationMethod();
        const MatrixType &DN_De = rGeometry.ShapeFunctionsLocalGradients(ThisIntegrationMethod)[rPointNumber];

        for (unsigned int i = 0; i < rGeometry.size(); i++)
        {
            const array_1d<double, 3> &coords = rGeometry[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
            for (unsigned int k = 0; k < rGeometry.WorkingSpaceDimension(); k++)
            {
                for (unsigned int m = 0; m < rGeometry.LocalSpaceDimension(); m++)
                {
                    rJ0(k, m) += coords[k] * DN_De(i, m);
                }
            }
        }

        KRATOS_CATCH("");
    }
    

    MatrixType CalculateShapeFunctionDerivatives(
    const int &rdimension, const double &rPointNumber, GeometryType& rGeometry)
    {
        KRATOS_TRY;

        const IntegrationMethod ThisIntegrationMethod = rGeometry.GetDefaultIntegrationMethod();
        uint dimension = rGeometry.WorkingSpaceDimension();
        uint number_of_nodes = rGeometry.size();
        MatrixType DN_DX(number_of_nodes, dimension);
        MatrixType J0(dimension, dimension);
        MatrixType invJ0(dimension, dimension);
        

        J0.clear();

        double detJ0;
        const MatrixType &DN_De = rGeometry.ShapeFunctionsLocalGradients(ThisIntegrationMethod)[rPointNumber];

        CalculateInitialJacobian(J0, rPointNumber, rGeometry);

        MathUtils<double>::InvertMatrix(J0, invJ0, detJ0);

        noalias(DN_DX) = prod(DN_De, invJ0);

        return DN_DX;

        KRATOS_CATCH("");
    }

    }  // namespace Move Mesh Utilities.
  
}  // namespace Kratos.