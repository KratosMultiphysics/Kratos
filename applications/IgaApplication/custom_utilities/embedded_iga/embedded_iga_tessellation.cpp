//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "embedded_iga_tessellation.h"

namespace Kratos
{
void EmbeddedIgaTessellation::CreateTessellation1D(
    const double mTessellationError,
    const BrepEdge& rCurveGeometry,
    std::vector<array_1d<double, 3>>& rPolygon)
{
    const int degree = rCurveGeometry.GetCurve3d()->Degree();
    const int number_cps = rCurveGeometry.GetCurve3d()->NbPoles();
    const bool is_rational = true;
    
    ANurbs::Pointer<ANurbs::CurveGeometry3D> curve_geometry = ANurbs::New<ANurbs::CurveGeometry3D>(degree, number_cps, is_rational);
    
    const std::vector<double> knot_vector = rCurveGeometry.GetCurve3d()->Knots();
    const unsigned int number_knots = rCurveGeometry.GetCurve3d()->NbKnots();

    for (unsigned int knot_i = 0; knot_i < number_knots; ++knot_i)
    {
        curve_geometry->SetKnot(knot_i, knot_vector[knot_i]);
    }
    for (unsigned int cp_i = 0; cp_i < number_cps; ++cp_i)
    {
        const auto node = rCurveGeometry.GetCurve3d()->GetNode(cp_i);
        curve_geometry->SetPole(cp_i, {node->X(), node->Y(), node->Z()});
        curve_geometry->SetWeight(cp_i, 
            rCurveGeometry.GetCurve3d()->GetNode(cp_i)->GetValue(NURBS_CONTROL_POINT_WEIGHT));  
    }

    // Create the three dimensional curve which is to be tessellated
    ANurbs::Curve3D curve(curve_geometry);

    // Perform the tessellation of the curve with flatness factor
    ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
    tessellation->Compute(curve, mTessellationError);

    rPolygon.resize(tessellation->NbPoints()); 

    for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)
    { 
        rPolygon[i][0] = tessellation->Point(i).X(); 
        rPolygon[i][1] = tessellation->Point(i).Y(); 
        rPolygon[i][2] = tessellation->Point(i).Z();
    }
}

void EmbeddedIgaTessellation::CreateTessellation2D(
    const double mTessellationError,
    const BrepFace& rFaceGeometry, 
    std::vector<std::vector<array_1d<double, 2>>>& rOuterPolygon,
    std::vector<std::vector<array_1d<double, 2>>>& rInnerPolygon)
{   
    const auto surface_geometry = rFaceGeometry.GetSurface();
    for (unsigned int b_loop_i = 0; b_loop_i < rFaceGeometry.GetBoundaryLoop().size(); ++b_loop_i)
    {
        std::vector<array_1d<double, 2>> loop;
        unsigned int point_index = 0;
        auto boundary_loop = rFaceGeometry.GetBoundaryLoop()[b_loop_i];

        for (unsigned int trim_i = 0; trim_i < boundary_loop.GetTrimmingCurves().size(); ++trim_i)
        {
            const auto trimming_curve = boundary_loop.GetTrimmingCurves()[trim_i];
            int trim_index = trimming_curve.GetTrimIndex();

            const auto curve_geometry = rFaceGeometry.GetTrimCurve(trim_index);

            const auto curve_on_surface = CurveOnSurface<3>(
                curve_geometry->CurveGeometry(),
                surface_geometry,
                curve_geometry->Domain());

            const auto tessellation = Kratos::make_shared<ANurbs::CurveTessellation<array_1d<double, 3>>>();

            tessellation->Compute(curve_on_surface, mTessellationError);

            // polygon vector needs to be resized to the current number of points + the new points
            // generated in the tessellation. However, the first point of one trimming curve
            // is the last point of the previous trimming curve. To account for these points
            // we subtract - 1.

            loop.resize(loop.size() + tessellation->NbPoints() - 1, ZeroVector(2));
                        
            for (unsigned int i = 0; i < tessellation->NbPoints() - 1; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    loop[point_index][j] = curve_geometry->CurveGeometry()
                                                        ->PointAt(tessellation->Parameter(i))[j];
                }
                point_index += 1;
            }
            
        }
        
        if (boundary_loop.IsOuterLoop())
        {
            rOuterPolygon.push_back(loop);
        }
        else 
        {
            rInnerPolygon.push_back(loop);
        }
    }
}

EmbeddedIgaTessellation::EmbeddedIgaTessellation()
{}

} // namespace Kratos.