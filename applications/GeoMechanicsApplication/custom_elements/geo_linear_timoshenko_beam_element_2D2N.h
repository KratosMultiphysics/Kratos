// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

#include "custom_elements/beam_elements/timoshenko_beam_element_2D2N.h"

namespace Kratos
{

/**
 * @class GeoLinearTimoshenkoBeamElement2D2N
 * @brief Geo-specific wrapper around LinearTimoshenkoBeamElement2D2N.
 *
 * This subclass overrides FinalizeSolutionStep to call
 * FinalizeMaterialResponsePK2 on each integration point's constitutive law,
 * enabling state-dependent laws (such as PiecewiseLinearMomentCapacityConstitutiveLaw)
 * to commit their internal state at the end of each converged time step.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearTimoshenkoBeamElement2D2N : public LinearTimoshenkoBeamElement2D2N
{
public:
    using BaseType = LinearTimoshenkoBeamElement2D2N;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoLinearTimoshenkoBeamElement2D2N);

    GeoLinearTimoshenkoBeamElement2D2N() = default;

    GeoLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    GeoLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoLinearTimoshenkoBeamElement2D2N>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoLinearTimoshenkoBeamElement2D2N>(NewId, pGeom, pProperties);
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

}; // class GeoLinearTimoshenkoBeamElement2D2N

} // namespace Kratos
