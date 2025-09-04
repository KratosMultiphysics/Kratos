//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "compute_beam_vectors_process.h"
#include "utilities/math_utils.h"
#include "includes/define.h"
#include "iga_application_variables.h"
#include "custom_elements/isogeometric_beam_element.h" 
namespace Kratos
{

//for normal beam    
ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(&(rModelPart.ElementsBegin()->GetGeometry().GetGeometryParent(0)))
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
    
    // Evaluate at the start of the beam (parameter = 0)
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    std::vector<array_1d<double, 3>> global_space_derivatives;
    
    mpParentCurve->GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);


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
    
    // Evaluate at the start of the beam (parameter = 0)
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    std::vector<array_1d<double, 3>> global_space_derivatives;
    mpParentCurve->GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);

    // Compute tangent vector (normalized)
    rT0 = global_space_derivatives[1];
    rT0 /= norm_2(rT0);

    // Compute normal vector (from surface)
    array_1d<double, 3> parameter_space_coordinates;
    const auto& nurbs_curve_2d = *(rCurve.pGetCurve());
    const auto& nurbs_surface_3d = *(rCurve.pGetSurface());
    nurbs_curve_2d.GlobalCoordinates(parameter_space_coordinates, local_coordinates);
    nurbs_surface_3d.GlobalSpaceDerivatives(global_space_derivatives, parameter_space_coordinates, 1);
    array_1d<double, 3> G1 = global_space_derivatives[1];
    array_1d<double, 3> G2 = global_space_derivatives[2];
    MathUtils<double>::CrossProduct(rN0, G1, G2); //surface normal = G1 × G2
    rN0 /= norm_2(rN0);
    
    MathUtils<double>::CrossProduct(V, rT0, rN0);
    V /= norm_2(V); 
}


void ComputeBeamVectorsProcess::ExecuteInitialize()
{
    const auto* curve_on_surface = dynamic_cast<const NurbsCurveOnSurfaceGeometry<3, PointerVector<Node>, PointerVector<Node>>*>(mpParentCurve);
    const auto* nurbs_curve = dynamic_cast<const NurbsCurveGeometry<3, PointerVector<Node>>*>(mpParentCurve);

    array_1d<double, 3> T;
    array_1d<double, 3> N;

    if (curve_on_surface != nullptr) {
        if (mpParentSurface != nullptr) {
            ComputeT0AndN0(*curve_on_surface, T, N);
        } else {
            ComputeT0AndN0(*nurbs_curve, T, N);
        }
        } else {
            KRATOS_ERROR << "Unknown geometry type in ComputeBeamVectorsProcess";
    }

    auto& properties = mrThisModelPart.GetProperties(0);
    properties.SetValue(T_0, T);
    properties.SetValue(N_0, N);
}


} // namespace Kratos