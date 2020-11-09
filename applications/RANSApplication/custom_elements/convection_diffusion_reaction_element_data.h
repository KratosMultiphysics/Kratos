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
#include "includes/node.h"
#include "includes/process_info.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class ConvectionDiffusionReactionElementData
{
public:
    using GeometryType = Geometry<Node<3>>;

    ConvectionDiffusionReactionElementData(const GeometryType& rGeometry)
    : mrGeometry(rGeometry)
    {
    }

    virtual void CalculateConstants(
        const ProcessInfo& rCurrentProcessInfo) = 0;

    virtual void CalculateGaussPointData(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const int Step) = 0;

    virtual array_1d<double, 3> CalculateEffectiveVelocity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateEffectiveKinematicViscosity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateReactionTerm(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual double CalculateSourceTerm(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const = 0;

    virtual void UpdateElementDataValueContainer(Element& rElement) const
    {
    }

    const GeometryType& GetGeometry() const
    {
        return mrGeometry;
    }

private:
    const GeometryType& mrGeometry;
};

///@}
} // namespace Kratos

#endif