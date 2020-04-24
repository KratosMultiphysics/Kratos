//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED

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
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class ScalarConvectionDiffusionReactionElementData
{
public:
    using GeometryType = Geometry<Node<3>>;

    ScalarConvectionDiffusionReactionElementData(const GeometryType& rGeometry)
        : mrGeometry(rGeometry)
    {
    }

    virtual void CalculateGaussPointData(const Vector& rShapeFunctions,
                                         const Matrix& rShapeFunctionDerivatives,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const int Step) = 0;

    virtual double CalculateEffectiveKinematicViscosity(
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateReactionTerm(const Vector& rShapeFunctions,
                                         const Matrix& rShapeFunctionDerivatives,
                                         const ProcessInfo& rCurrentProcessInfo) const = 0;

    virtual double CalculateSourceTerm(const Vector& rShapeFunctions,
                                       const Matrix& rShapeFunctionDerivatives,
                                       const ProcessInfo& rCurrentProcessInfo) const = 0;

    const GeometryType& GetGeometry() const
    {
        return mrGeometry;
    }

private:
    const GeometryType& mrGeometry;
};
} // namespace Kratos

#endif