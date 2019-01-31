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
    void EmbeddedIgaTessellation::CreateTessellationCurveOnSurface(std::vector<array_1d<double, 2> >& rPolygon)
    {
        unsigned int point_id = 0; 
        for (unsigned int brep_i = 0; brep_i < mBrepModelVector.size(); ++brep_i)
        {
            for (unsigned int face_i = 0; face_i < mBrepModelVector[brep_i].GetFaceVector().size(); ++face_i)
            {
                auto surface_geometry = mBrepModelVector[brep_i].GetFaceVector()[face_i].GetSurface();

                for (unsigned int b_loop_i = 0; b_loop_i < mBrepModelVector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop().size(); ++b_loop_i)
                {   
                    for (unsigned int trim_i = 0; trim_i < mBrepModelVector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop()[b_loop_i].GetTrimmingCurves().size(); ++trim_i)
                    {   
                        auto trimming_curves = mBrepModelVector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop()[b_loop_i].GetTrimmingCurves();
                        int trim_index = trimming_curves[trim_i].GetTrimIndex(); 

                        auto curve_geometry = mBrepModelVector[brep_i].GetFaceVector()[face_i].GetTrimCurve(
                            trim_index); 

                        auto curve_on_surface = CurveOnSurface<3>(
                            curve_geometry->CurveGeometry(), 
                            surface_geometry, 
                            curve_geometry->Domain());

                        auto tessellation = Kratos::make_shared<ANurbs::CurveTessellation<array_1d<double, 3>>>();
                        
                        tessellation->Compute(curve_on_surface, 1e-4);
            
                        rPolygon.resize(point_id + tessellation->NbPoints() - 1); 
        
                        for (unsigned int i = 0; i < tessellation->NbPoints() -1; ++i)
                        {
                            double curve_para = tessellation->Parameter(i); 
                            for (unsigned int j = 0; j < 2; ++j)    rPolygon[point_id][j] = curve_geometry->CurveGeometry()->PointAt(curve_para)[j];      
                            point_id += 1; 
                        }
                    }
                }
            }
        }
    }

    void EmbeddedIgaTessellation::CreateTessellationCurve(std::vector<array_1d<double, 3> >& rPolygon)
    {
        unsigned int point_id = 0; 
        for (unsigned int brep_i = 0; brep_i < mBrepModelVector.size(); ++brep_i)
        {   
            for (unsigned int edge_i = 0; edge_i < mBrepModelVector[brep_i].GetEdgeVector().size(); ++edge_i)
            {
                const int degree = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Degree();
                const int number_cps = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbPoles();
                const bool is_rational = true;

                ANurbs::Pointer<ANurbs::CurveGeometry3D> geometry = ANurbs::New<ANurbs::CurveGeometry3D>(degree, number_cps, is_rational);
                
                const std::vector<double> knot_vector = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Knots();
                const unsigned int number_knots = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbKnots();

                for (unsigned int knot_i = 0; knot_i < number_knots; ++knot_i)
                {
                    geometry->SetKnot(knot_i, knot_vector[knot_i]);
                }
                for (unsigned int cp_i = 0; cp_i < number_cps; ++cp_i)
                {
                    const auto node = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(cp_i);
                    geometry->SetPole(cp_i, {node->X(), node->Y(), node->Z()});
                    geometry->SetWeight(cp_i, mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(cp_i)->GetValue(NURBS_CONTROL_POINT_WEIGHT));  
                }
    
                // Create the three dimensional curve which is to be tessellated
                ANurbs::Curve3D curve(geometry);

                // Perform the tessellation of the curve with flatness factor
                ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
                
                tessellation->Compute(curve, 1e-4);

                rPolygon.resize(point_id + tessellation->NbPoints() - 1); 
            
                int trim_index = mBrepModelVector[brep_i].GetEdgeVector()[edge_i].GetBrepEdgeTopologyVector()[0].trim_index; 
                for (unsigned int i = 0; i < tessellation->NbPoints() - 1; ++i)
                {
                    double para = tessellation->Parameter(i);

                    auto point_UV = mBrepModelVector[brep_i].GetFaceVector()[0].GetTrimCurve(trim_index)->CurveGeometry()->PointAt(para); 

                    KRATOS_WATCH(point_UV)

                    rPolygon[point_id][0] = point_UV[0]; 
                    rPolygon[point_id][1] = point_UV[1]; 
                    rPolygon[point_id][2] = 0; 
                    point_id += 1; 
                }


                // for (unsigned int i = 0; i < tessellation->NbPoints() -1; ++i)
                // {
                //     rPolygon[point_id][0] = tessellation->Point(i).X(); 
                //     rPolygon[point_id][1] = tessellation->Point(i).Y(); 
                //     rPolygon[point_id][2] = tessellation->Point(i).Z();
                //     point_id += 1; 
                // }
            }
        }
    }

    void EmbeddedIgaTessellation::CreateTessellationParameterCurve(std::vector<array_1d<double, 3> >& rPolygon)
    {
        unsigned int point_id = 0; 
        for (unsigned int brep_i; brep_i < mBrepModelVector.size(); ++brep_i)
        {
            for (unsigned int face_i = 0; face_i < mBrepModelVector[brep_i].GetFaceVector().size(); ++face_i)
            {   
                auto boundary_loop = mBrepModelVector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop();
        
                for (unsigned int b_loop_i = 0; b_loop_i < boundary_loop.size(); ++b_loop_i)
                {
                    for (unsigned int t_curve_i = 0; t_curve_i < boundary_loop[b_loop_i].GetTrimmingCurves().size(); ++t_curve_i)
                    {
                        const unsigned int degree = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->Degree();
                        const unsigned int number_cps = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->NbPoles();
                        const bool is_rational = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->IsRational();
                        const unsigned int number_knots = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->NbKnots();
                        const std::vector<double> knot_vector = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Knots();
                        
                        ANurbs::Pointer<ANurbs::CurveGeometry2D> geometry = ANurbs::New<ANurbs::CurveGeometry2D>(degree, number_cps, is_rational);

                        for (unsigned int knot_i = 0; knot_i < number_knots; ++knot_i)
                        {
                            geometry->SetKnot(knot_i, knot_vector[knot_i]);
                        }
                        
                        for (unsigned int cp_i = 0; cp_i < number_cps; ++cp_i)
                        {
                            geometry->SetPole(cp_i, {boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Poles()[cp_i][0], 
                                                     boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Poles()[cp_i][1]});
                            if (is_rational == true)
                            {   
                                geometry->SetWeight(cp_i, boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Weight(cp_i)); 
                            }
                        }

                        // Create the three dimensional curve which is to be tessellated
                        ANurbs::Curve2D curve(geometry);

                        // Perform the tessellation of the curve with flatness factor 
                        auto tessellation = ANurbs::New<ANurbs::CurveTessellation2D>();
                        tessellation->Compute(curve, 1e-2);
                        
                        rPolygon.resize(point_id + tessellation->NbPoints() - 1); 
                        
                        for (unsigned int i = 0; i < tessellation->NbPoints() -1; ++i)
                        {
                            rPolygon[point_id][0] = tessellation->Point(i).X(); 
                            rPolygon[point_id][1] = tessellation->Point(i).Y(); 
                            // rPolygon[point_id][2] = tessellation->Point(i).Z();
                            rPolygon[point_id][2] = 0;
                            point_id += 1; 
                        }
                    }
                }
            }
        }
    }
    EmbeddedIgaTessellation::EmbeddedIgaTessellation(std::vector<BrepModel>&  rBrepModelVector)
        : mBrepModelVector(rBrepModelVector)
    {
    }
} // namespace Kratos.