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

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class ConvectionDiffusionReactionElementData
{
public:
    using GeometryType = Geometry<Node<3>>;

    ConvectionDiffusionReactionElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
        : mrGeometry(rGeometry),
          mrProperties(rProperties),
          mrConstitutiveLaw(*(rGeometry.GetValue(CONSTITUTIVE_LAW)))
    {
        mConstitutiveLawParameters =
            ConstitutiveLaw::Parameters(rGeometry, rProperties, rProcessInfo);
    }

    virtual void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
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

///@}
} // namespace Kratos

#endif