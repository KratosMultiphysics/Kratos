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
void EmbeddedIgaTessellation::CreateTessellation(
    const BrepFace& rFaceGeometry, 
    std::vector<std::vector<array_1d<double, 2>>>& rOuterPolygon,
    std::vector<std::vector<array_1d<double, 2>>>& rInnerPolygon)
{   
    std::vector<array_1d<double, 2>> loop;
    const auto surface_geometry = rFaceGeometry.GetSurface();
    for (unsigned int b_loop_i = 0; b_loop_i < rFaceGeometry.GetBoundaryLoop().size(); ++b_loop_i)
    {
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

            tessellation->Compute(curve_on_surface, 1e-3);


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