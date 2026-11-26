//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti,
//                   Ricky Aristio, Alicia Knauer
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "compute_beam_vectors_process.h"
#include "utilities/math_utils.h"
#include "includes/define.h"
#include "iga_application_variables.h"
#include "custom_elements/non_linear_bernoulli_beam_element_3D.h"
#include "custom_elements/bernoulli_beam_base_element_3D.h"

namespace Kratos
{

//for normal beam
ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(nullptr)
{
}

//for normal beam
ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart,
    const NurbsCurveGeometry<3, PointerVector<Node>>& rParentCurve)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(&rParentCurve)
{
}

//for embedded beam
ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart,
    const NurbsCurveOnSurfaceGeometry<3, PointerVector<Node>, PointerVector<Node>>& rParentCurve)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(&rParentCurve)
{
}

void ComputeBeamVectorsProcess::ComputeT0AndN0(
    const NurbsCurveGeometry<3, PointerVector<Node>>& rCurve,
    array_1d<double, 3>& rT0,
    array_1d<double, 3>& rN0)
{
    array_1d<double, 3> V = ZeroVector(3);
    array_1d<double, 3> up = ZeroVector(3);

    // Evaluate at the start of the beam
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    local_coordinates[0] = rCurve.DomainInterval().GetT0();
    std::vector<array_1d<double, 3>> global_space_derivatives;

    rCurve.GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);


    // Tangent vector (normalized)
    rT0 = global_space_derivatives[1];
    rT0 /= norm_2(rT0);

    up[2] = 1.0; // Z-up
    if (std::abs(inner_prod(rT0, up)) > 0.99) {
        up[2] = 0.0;
        up[1] = 1.0; // switch to Y-up
    }

    MathUtils<double>::CrossProduct(V, rT0, up);
    V /= norm_2(V);
    MathUtils<double>::CrossProduct(rN0, V, rT0);
    rN0 /= norm_2(rN0);
}

void ComputeBeamVectorsProcess::ComputeT0AndN0(
    const NurbsCurveOnSurfaceGeometry<3, PointerVector<Node>, PointerVector<Node>>& rCurve,
    array_1d<double, 3>& rT0,
    array_1d<double, 3>& rN0)
{
    array_1d<double, 3> V = ZeroVector(3);
    array_1d<double, 3> up = ZeroVector(3);

    // Evaluate at the start of the beam
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    local_coordinates[0] = rCurve.DomainInterval().GetT0();
    std::vector<array_1d<double, 3>> global_space_derivatives;
    rCurve.GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);

    // Compute tangent vector (normalized)
    rT0 = global_space_derivatives[1];
    rT0 /= norm_2(rT0);

    // Compute normal vector (from surface)
    array_1d<double, 3> parameter_space_coordinates;
    const auto& nurbs_curve_2d = *(rCurve.pGetCurve());
    const auto p_surface_geometry = rCurve.pGetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX);
    const auto& nurbs_surface_3d = static_cast<const NurbsSurfaceGeometry<3, Geometry<Node>::PointsArrayType>&>(*p_surface_geometry);
    nurbs_curve_2d.GlobalCoordinates(parameter_space_coordinates, local_coordinates);
    nurbs_surface_3d.GlobalSpaceDerivatives(global_space_derivatives, parameter_space_coordinates, 1);
    array_1d<double, 3> G1 = global_space_derivatives[1];
    array_1d<double, 3> G2 = global_space_derivatives[2];
    MathUtils<double>::CrossProduct(rN0, G1, G2); //surface normal = G1 × G2
    rN0 /= norm_2(rN0);

    MathUtils<double>::CrossProduct(V, rT0, rN0);
    V /= norm_2(V);
}


void ComputeBeamVectorsProcess::ComputeT0AndN0(
    const Geometry<Node>& rParentCurve,
    array_1d<double, 3>& rT0,
    array_1d<double, 3>& rN0)
{
    const auto* curve_on_surface = dynamic_cast<const NurbsCurveOnSurfaceGeometry<3, PointerVector<Node>, PointerVector<Node>>*>(&rParentCurve);
    auto* nurbs_curve = dynamic_cast<const NurbsCurveGeometry<3, PointerVector<Node>>*>(&rParentCurve);
    const auto* brep_curve = dynamic_cast<const BrepCurve<PointerVector<Node>, PointerVector<Point>>*>(&rParentCurve);

    if (brep_curve != nullptr) {
        nurbs_curve = dynamic_cast<const NurbsCurveGeometry<3, PointerVector<Node>>*>(
            &brep_curve->GetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX));
    }

    if (curve_on_surface != nullptr) {
        ComputeT0AndN0(*curve_on_surface, rT0, rN0);
    } else if (nurbs_curve != nullptr) {
        ComputeT0AndN0(*nurbs_curve, rT0, rN0);
    } else {
        KRATOS_ERROR << "Unknown geometry type in ComputeBeamVectorsProcess";
    }
}

void ComputeBeamVectorsProcess::ExecuteInitialize()
{
    // Each beam element gets the initial frame of its own parent curve, so that several beams can share a model part.
    std::unordered_map<const Geometry<Node>*, std::pair<array_1d<double, 3>, array_1d<double, 3>>> frames;

    for (auto& r_element : mrThisModelPart.Elements()) {
        const Geometry<Node>* p_parent_curve = &r_element.GetGeometry().GetGeometryParent(0);

        if (mpParentCurve != nullptr && p_parent_curve != mpParentCurve) {
            continue;
        }

        auto it_frame = frames.find(p_parent_curve);
        if (it_frame == frames.end()) {
            array_1d<double, 3> T;
            array_1d<double, 3> N;
            ComputeT0AndN0(*p_parent_curve, T, N);
            it_frame = frames.emplace(p_parent_curve, std::make_pair(T, N)).first;
        }

        r_element.SetValue(T_0, Vector(it_frame->second.first));
        r_element.SetValue(N_0, Vector(it_frame->second.second));
    }
}


} // namespace Kratos
