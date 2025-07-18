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

namespace Kratos
{

ComputeBeamVectorsProcess::ComputeBeamVectorsProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
    , mrThisModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
{
}

void ComputeBeamVectorsProcess::ExecuteInitialize()
{
    for (auto& r_element : mrThisModelPart.Elements()) {
        auto& r_geometry = r_element.GetGeometry();
        
        // Check if this is a NURBS curve geometry (beam element)
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Nurbs_Curve) {
            
            // Get the NURBS curve geometry - handle both direct curves and quadrature point geometries
            const NurbsCurveGeometry<3, PointerVector<Node>>* p_curve = nullptr;
            
            // Try to get parent geometry first (for quadrature point geometries)
            try {
                auto& parent_geometry = r_geometry.GetGeometryParent(0);
                p_curve = dynamic_cast<const NurbsCurveGeometry<3, PointerVector<Node>>*>(&parent_geometry);
            } catch (...) {
                // Direct NURBS curve geometry
                p_curve = dynamic_cast<const NurbsCurveGeometry<3, PointerVector<Node>>*>(&r_geometry);
            }
            
            if (p_curve != nullptr) {
                // Compute T0 and N0 vectors
                array_1d<double, 3> t0_vector, n0_vector;
                this->ComputeT0AndN0(*p_curve, t0_vector, n0_vector);
                
                // Convert to Kratos Vector format and set as properties
                Vector t0_kratos(3), n0_kratos(3);
                for (int i = 0; i < 3; ++i) {
                    t0_kratos[i] = t0_vector[i];
                    n0_kratos[i] = n0_vector[i];
                }
                
                // Set the computed vectors as element properties
                auto& r_properties = r_element.GetProperties();
                r_element.SetValue(T_0, t0_kratos);
                r_element.SetValue(N_0, n0_kratos);
            }
        }
    }
}

void ComputeBeamVectorsProcess::ComputeT0AndN0(
    const NurbsCurveGeometry<3, PointerVector<Node>>& rCurve,
    array_1d<double, 3>& rT0,
    array_1d<double, 3>& rN0)
{
    // Set up local coordinates [0, 0, 0] (parameter space)
    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 0.0;
    local_coordinates[1] = 0.0;
    local_coordinates[2] = 0.0;
    
    // Prepare container for derivatives
    std::vector<array_1d<double, 3>> global_space_derivatives;
    
    // Compute global space derivatives up to order 1
    rCurve.GlobalSpaceDerivatives(global_space_derivatives, local_coordinates, 1);
    
    // Extract the first derivative (tangent vector) - it's at index 1
    array_1d<double, 3> t0_unnormalized = global_space_derivatives[1];
    
    // Normalize T0
    double t0_norm = norm_2(t0_unnormalized);
    rT0 = t0_unnormalized / t0_norm;
    
    // Compute N0 = cross([0, -1, 0], t0)
    array_1d<double, 3> y_axis_negative;
    y_axis_negative[0] = 0.0;
    y_axis_negative[1] = -1.0;
    y_axis_negative[2] = 0.0;
    
    MathUtils<double>::CrossProduct(rN0, y_axis_negative, rT0);
}



} // namespace Kratos