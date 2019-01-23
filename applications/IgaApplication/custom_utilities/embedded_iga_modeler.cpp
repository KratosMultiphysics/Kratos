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
#include "embedded_iga_modeler.h"


namespace Kratos
{

    void EmbeddedIgaModeler::CreateTessellationCurve(ANurbs::Pointer<ANurbs::CurveTessellation3D>& rTessellation)
    {
        for (unsigned int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
        {
            for (unsigned int edge_i = 0; edge_i < m_brep_model_vector[brep_i].GetEdgeVector().size(); ++edge_i)
            {
                const int degree = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Degree();
                const int number_cps = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbPoles();
                const bool is_rational = true;

                ANurbs::Pointer<ANurbs::CurveGeometry3D> geometry = ANurbs::New<ANurbs::CurveGeometry3D>(degree, number_cps, is_rational);
                
                const std::vector<double> knot_vector = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Knots();
                const int number_knots = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbKnots();

                for (unsigned int knot_i = 0; knot_i < number_knots; ++knot_i)
                {
                    geometry->SetKnot(knot_i, knot_vector[knot_i]);
                }

                for (unsigned int cp_i = 0; cp_i < number_cps; ++cp_i)
                {
                    const auto node = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(cp_i);
                    const double weight = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(cp_i)->GetValue(NURBS_CONTROL_POINT_WEIGHT); 
    
                    geometry->SetPole(cp_i, {node->X(), node->Y(), node->Z()});
                    geometry->SetWeight(cp_i, weight);  
                }
    
                // Create the three dimensional curve which is to be tessellated
                ANurbs::Curve3D curve(geometry);
                // Perform the tessellation of the curve with flatness factor 1e-2
                rTessellation->Compute(curve, 1e-5);
            }
        }
    }

    void EmbeddedIgaModeler::CreateElements2D(ModelPart& rSkinModelPart)
    {
        ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
        CreateTessellationCurve(tessellation);

        // Create Nodes in Skin Modelpart from the Points generated in the tessellation of the curve
        for (unsigned int point_i = 0; point_i < tessellation->NbPoints(); ++point_i)
        {
            rSkinModelPart.CreateNewNode(point_i, tessellation->Point(point_i).X(), tessellation->Point(point_i).Y(), tessellation->Point(point_i).Z());
        }

        // skin_model_part property Container (Container is empty, but for the skin no properties are needed)
        Properties::Pointer prop = rSkinModelPart.pGetProperties(0);

        unsigned int node_id = 0;
        // Create Elements in skin_model_part
        for (unsigned int element_i = 0; element_i < tessellation->NbPoints() - 1; ++element_i)
        {
            rSkinModelPart.CreateNewElement("Element2D2N", element_i, {{node_id, node_id + 1}}, prop);
            node_id += 1;
        }
    }

    void EmbeddedIgaModeler::CreateTessellationParameterCurve(std::vector<array_1d<double, 3> >& rPolygon)
    {
        unsigned int point_id = 0; 
        for (unsigned int brep_i; brep_i < m_brep_model_vector.size(); ++brep_i)
        {
            for (unsigned int face_i = 0; face_i < m_brep_model_vector[brep_i].GetFaceVector().size(); ++face_i)
            {   
                auto boundary_loop = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop();
        
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
                        // Perform the tessellation of the curve with flatness factor 1e-2
                        ANurbs::Pointer<ANurbs::CurveTessellation2D> tessellation = ANurbs::New<ANurbs::CurveTessellation2D>();
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
        
        rPolygon.resize(point_id + 1);
        rPolygon[point_id][0] = rPolygon[0][0];
        rPolygon[point_id][1] = rPolygon[0][1];
        rPolygon[point_id][2] = rPolygon[0][2];
    }

    std::vector<Matrix> EmbeddedIgaModeler::Triangulate()
    {
        std::vector<array_1d<double,3>> polygon;

        CreateTessellationParameterCurve(polygon);
        
        std::vector<Matrix> triangles;
        EmbeddedIgaTriangulation embedded_triangulation; 
        embedded_triangulation.CreateTrianglesEmpire(polygon, triangles); 
        
        return triangles;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintCurveTessellationPoints()
    {
        ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
        CreateTessellationCurve(tessellation);
        std::vector<std::vector<double>> coords(tessellation->NbPoints(), std::vector<double>(2,0));

        for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)
        {    
            coords[i][0] = tessellation->Point(i).X();
            coords[i][1] = tessellation->Point(i).Y();
        }
        return coords; 
    }

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintTriangulationPoints()
    {
        std::vector<Matrix> triangles = Triangulate();
        std::vector<std::vector<double>> coords(triangles.size()*3, std::vector<double>(2.0)); 
        
        unsigned int index = 0; 

        for (unsigned int i = 0; i < triangles.size(); ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                coords[index][0] = triangles[i](j,0); 
                coords[index][1] = triangles[i](j,1); 
                index += 1; 
            }
        }
        return coords; 
    }

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintParameterCurveTessellationPoints()
    {
        std::vector<array_1d<double,3>> polygon; 
        CreateTessellationParameterCurve(polygon);
        
        std::vector<std::vector<double>> coords(polygon.size(), std::vector<double>(2,0));

        for (unsigned int i = 0; i < polygon.size(); ++i)
        {
            coords[i][0] = polygon[i][0]; 
            coords[i][1] = polygon[i][1]; 
        }
        return coords;
    }

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintGaussPoints()
    {
        std::vector<Matrix> triangles = Triangulate(); 
        EmbeddedIgaErrorEstimation error_estimator(triangles); 
        std::vector<array_1d<double, 2> > gp_pos; 
        error_estimator.InsertGaussPoints(gp_pos);

        std::vector<std::vector<double>> gp_coords(gp_pos.size(), std::vector<double>(2,0)); 
        
        for (unsigned int i = 0; i < gp_pos.size(); ++i)
        {
            gp_coords[i][0] = gp_pos[i][0]; 
            gp_coords[i][1] = gp_pos[i][1]; 
        }
        return gp_coords; 
    }

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintMappedGaussPoints()
    {
        std::vector<Matrix> triangles = Triangulate(); 
        EmbeddedIgaErrorEstimation error_estimator(triangles); 
        std::vector<array_1d<double, 2> > gp_pos; 

        error_estimator.InsertGaussPoints(gp_pos); 
        

        auto geometry = m_brep_model_vector[0].GetFaceVector()[0].GetSurface(); 

        std::vector<std::vector<double>> mapped_gp_coords (gp_pos.size(), std::vector<double>(3,0)); 

        for (unsigned int i = 0; i < gp_pos.size(); ++i)
        {    
            auto point = geometry->PointAt(gp_pos[i][0],gp_pos[i][1]); 
            
            mapped_gp_coords[i][0] = point[0]; 
            mapped_gp_coords[i][1] = point[1]; 
            mapped_gp_coords[i][2] = point[2]; 
        }
        return mapped_gp_coords; 
    }


    void EmbeddedIgaModeler::TestTriangle()
    {
        std::vector<array_1d<double,3>> polygon;

        CreateTessellationParameterCurve(polygon);

    }

    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
        : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
        m_model_part(rModelPart)
    {
    }
} // namespace Kratos.