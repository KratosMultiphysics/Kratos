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
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
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
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo,
        ConstitutiveLaw& rConstitutiveLaw)
        : mrGeometry(rGeometry),
          mrProperties(rProperties),
          mrConstitutiveLaw(rConstitutiveLaw)
    {
        mConstitutiveLawParameters =
            ConstitutiveLaw::Parameters(rGeometry, rProperties, rProcessInfo);
    }

    ConstitutiveLaw::Parameters& GetConstitutiveLawParameters()
    {
        return mConstitutiveLawParameters;
    }

    ConstitutiveLaw& GetConstitutiveLaw()
    {
        return mrConstitutiveLaw;
    }

    const GeometryType& GetGeometry() const
    {
        return mrGeometry;
    }

    const Properties& GetProperties() const
    {
        return mrProperties;
    }

private:
    const GeometryType& mrGeometry;
    const Properties& mrProperties;
    ConstitutiveLaw& mrConstitutiveLaw;
    ConstitutiveLaw::Parameters mConstitutiveLawParameters;
};
} // namespace Kratos

#endif