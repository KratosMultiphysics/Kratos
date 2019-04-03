//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                       license: license.txt
//
//  License:          BSD License
//  Main authors:     Riccardo Rossi
//


#if !defined(KRATOS_INTEGRATION_UTILITIES_INCLUDED )
#define  KRATOS_INTEGRATION_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"

namespace Kratos
{
class IntegrationUtilities
{
public:
    static Geometry<Node<3> >::IntegrationMethod GetIntegrationMethodForExactMassMatrixEvaluation(Geometry<Node<3>> const& geom)
    {
        Geometry<Node<3> >::IntegrationMethod integration_method = geom.GetDefaultIntegrationMethod();
        if(integration_method == GeometryData::GI_GAUSS_1)
            integration_method = GeometryData::GI_GAUSS_2;
        else if(integration_method == GeometryData::GI_GAUSS_2)
            integration_method = GeometryData::GI_GAUSS_3;
        else if(integration_method == GeometryData::GI_GAUSS_3)
            integration_method = GeometryData::GI_GAUSS_4;
        return integration_method;
    }
};

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_UTILITIES_INCLUDED  defined


