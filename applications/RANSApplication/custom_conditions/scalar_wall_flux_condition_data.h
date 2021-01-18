//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_SCALAR_WALL_FLUX_CONDITION_DATA_H_INCLUDED)
#define KRATOS_SCALAR_WALL_FLUX_CONDITION_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class ScalarWallFluxConditionData
{
public:
    using GeometryType = Geometry<Node<3>>;

    ScalarWallFluxConditionData(
        const GeometryType& rGeometry)
    : mrGeometry(rGeometry)
    {
    }

    virtual void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo) = 0;

    virtual bool IsWallFluxComputable() const = 0;

    virtual double CalculateWallFlux(
        const Vector& rShapeFunctions) const = 0;

    const GeometryType& GetGeometry() const
    {
        return mrGeometry;
    }

private:
    const GeometryType& mrGeometry;
};
} // namespace Kratos

#endif