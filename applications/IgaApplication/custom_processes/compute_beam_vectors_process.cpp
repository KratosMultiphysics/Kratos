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
    , mrParentCurve(rParentCurve)
{
}

void ComputeBeamVectorsProcess::ExecuteInitialize()
{
    // Use the provided parent curve
    array_1d<double, 3> t0, n0;
    ComputeT0AndN0(mrParentCurve, t0, n0);

    // Set T_0 and N_0 on the properties (elements access via GetProperties()[T_0/N_0])
    auto& properties = mrThisModelPart.GetProperties(0);
    properties.SetValue(T_0, t0);
    properties.SetValue(N_0, n0);
}

void ComputeBeamVectorsProcess::ComputeT0AndN0(
    const NurbsCurveGeometry<3, PointerVector<Node>>& rCurve,
    array_1d<double, 3>& rT0,
    array_1d<double, 3>& rN0)
{
    // Evaluate at the start of the beam (parameter = 0)
    array_1d<double, 3> local_coordinates = ZeroVector(3);
    std::vector<array_1d<double, 3>> global_space_derivatives;
    rCurve.GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);

    // Tangent vector (normalized)
    rT0 = global_space_derivatives[1];
    rT0 /= norm_2(rT0);

    // Normal vector (cross product of -Y axis and tangent)
    array_1d<double, 3> y_axis_negative = ZeroVector(3);
    y_axis_negative[1] = -1.0;

    MathUtils<double>::CrossProduct(rN0, y_axis_negative, rT0);
}

} // namespace Kratos