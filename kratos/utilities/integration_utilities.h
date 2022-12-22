//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
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
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_4)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_5;
        return integration_method;
    }
};

}  // namespace Kratos.
