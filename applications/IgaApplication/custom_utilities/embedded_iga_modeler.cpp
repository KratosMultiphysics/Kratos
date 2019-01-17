//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Dagmawi Bekel
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
            rTessellation->Compute(curve, 1e-2);
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
    // needed later
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
                        const auto node_x = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Poles()[cp_i][0];
                        const auto node_y = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Poles()[cp_i][1];
                        const double weight = boundary_loop[b_loop_i].GetTrimmingCurves()[t_curve_i].GetCurve2D()->CurveGeometry()->Weight(cp_i); 

                        geometry->SetPole(cp_i, {node_x, node_y});
                        if (is_rational == true)   geometry->SetWeight(cp_i, weight); 
                    }

                    // Create the three dimensional curve which is to be tessellated
                    ANurbs::Curve2D curve(geometry);
                    // Perform the tessellation of the curve with flatness factor 1e-2
                    ANurbs::Pointer<ANurbs::CurveTessellation2D> tessellation = ANurbs::New<ANurbs::CurveTessellation2D>();
                    tessellation->Compute(curve, 1e-2);
                    
                    rPolygon.resize(point_id + tessellation->NbPoints()); 
                    
                    for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)
                    {
                        rPolygon[point_id][0] = tessellation->Point(i).X(); 
                        rPolygon[point_id][1] = tessellation->Point(i).Y(); 
                        // rPolygon[point_id][2] = tessellation->Point(i).Z(); // // 
                        rPolygon[point_id][2] = 0; 
                        point_id += 1; 
                    }
                }
            }
        }
    }
}


void EmbeddedIgaModeler::CreateElements3D()
{
    std::vector<array_1d<double,3>> polygon;
  

    CreateTessellationParameterCurve(polygon); 
    
    // Triangulation

    const unsigned int number_points = polygon.size(); 
    KRATOS_WATCH(number_points)
    
    
}



std::vector<double> EmbeddedIgaModeler::PrintNodesX()
{
    ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
    CreateTessellationCurve(tessellation);
    std::vector<double> x(tessellation->NbPoints());
    for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)     x[i] = tessellation->Point(i).X();
    return x;
}

std::vector<double> EmbeddedIgaModeler::PrintNodesY()
{
    ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>();
    CreateTessellationCurve(tessellation);
    std::vector<double> y(tessellation->NbPoints());
    for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)    y[i] = tessellation->Point(i).Y();
    return y;
}

std::vector<double> EmbeddedIgaModeler::PrintNodesX3D()
{
    // std::vector<double> x;
    // std::vector<double> y; 
    // CreateTessellationParameterCurve(x,y);
    
    // return x;
}

std::vector<double> EmbeddedIgaModeler::PrintNodesY3D()
{
    // std::vector<double> x;
    // std::vector<double> y; 
    // CreateTessellationParameterCurve(x,y);

    // KRATOS_WATCH(x)
    // KRATOS_WATCH(y)
    
    // return y;
}





EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
      m_model_part(rModelPart)
{
}
} // namespace Kratos.
