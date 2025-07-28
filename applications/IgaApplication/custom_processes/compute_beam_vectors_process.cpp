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


ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart,
    const NurbsCurveGeometry<3, PointerVector<Node>>& rParentCurve)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(&rParentCurve)
{
}

ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    ModelPart& rModelPart,
    const NurbsCurveOnSurfaceGeometry<3, PointerVector<Node>, PointerVector<Node>>& rParentCurve)
    : Process()
    , mrThisModelPart(rModelPart)
    , mpParentCurve(&rParentCurve)
{
}

void ComputeBeamVectorsProcess::ExecuteInitialize()
{
    array_1d<double, 3> T;
    array_1d<double, 3> V = ZeroVector(3);
    array_1d<double, 3> N = ZeroVector(3);
    array_1d<double, 3> up = ZeroVector(3);
    
    // Evaluate at the start of the beam (parameter = 0)
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    std::vector<array_1d<double, 3>> global_space_derivatives;
    mpParentCurve->GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);

    // Tangent vector (normalized)
    T = global_space_derivatives[1];
    T /= norm_2(T);

    up[2] = 1.0; // Z-up
    if (std::abs(inner_prod(T, up)) > 0.99) {
        up[2] = 0.0;
        up[1] = 1.0; // switch to Y-up
    }

    MathUtils<double>::CrossProduct(V, T, up);
    V /= norm_2(V); // normalize binormal
    MathUtils<double>::CrossProduct(N, V, T);
    N /= norm_2(N); 

    // Set T_0 and N_0 on the properties (elements access via GetProperties()[T_0/N_0])
    auto& properties = mrThisModelPart.GetProperties(0);
    properties.SetValue(T_0, T);
    properties.SetValue(N_0, N);

}


} // namespace Kratos