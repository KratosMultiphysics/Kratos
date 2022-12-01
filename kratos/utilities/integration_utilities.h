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

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"

namespace Kratos
{
class IntegrationUtilities
{
public:
    template<class TPointType>
    static GeometryData::IntegrationMethod GetIntegrationMethodForExactMassMatrixEvaluation(Geometry<TPointType> const& rGeometry)
    {
        GeometryData::IntegrationMethod integration_method = rGeometry.GetDefaultIntegrationMethod();
        if (integration_method == GeometryData::IntegrationMethod::GI_GAUSS_1)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_2)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_3;
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_3)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
        return integration_method;
    }
};

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_UTILITIES_INCLUDED  defined
