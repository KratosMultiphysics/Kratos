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

void EmbeddedIgaModeler::CreateTessellation2D(ANurbs::Pointer<ANurbs::CurveTessellation3D>& tessellation)
{   
    for (unsigned int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
    {
        for (unsigned int edge_i = 0; edge_i < m_brep_model_vector[brep_i].GetEdgeVector().size(); ++edge_i)
        {   
            const int degree = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Degree();
            const int number_cps = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbPoles();
            const bool isRational = false; 

            ANurbs::Pointer<ANurbs::CurveGeometry3D> geometry = ANurbs::New<ANurbs::CurveGeometry3D>(degree, number_cps, isRational);
        
            const std::vector<double> knot_vector = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Knots(); 
            const int number_knots = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbKnots(); 
            
            for (unsigned int knot_i = 0; knot_i < number_knots; ++knot_i)  
            {
                geometry->SetKnot(knot_i,knot_vector[knot_i]); 
            }
    
            for (unsigned int cp_i = 0; cp_i < number_cps; ++cp_i)
            {
                const auto node = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(cp_i); 
                geometry->SetPole(cp_i,{node->X(), node->Y(), node->Z()});
            }

            // Create the three dimensional curve which is to be tessellated
            ANurbs::Curve3D curve(geometry); 
            // Perform the tessellation of the curve with flatness factor 1e-2
            tessellation->Compute(curve, 1e-3); 
            
        }
    }
}

void EmbeddedIgaModeler::CreateElements2D(ModelPart& rSkinModelPart)
{
    
    ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>(); 
    CreateTessellation2D(tessellation);  
    
    // Create Nodes in Skin Modelpart from the Points generated in the tessellation of the curve    
    for (unsigned int point_i = 0; point_i < tessellation->NbPoints(); ++point_i)    
    {   
        rSkinModelPart.CreateNewNode(point_i, tessellation->Point(point_i).X(), tessellation->Point(point_i).Y(),tessellation->Point(point_i).Z());
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

// std::vector<double> EmbeddedIgaModeler::PrintNodesX()
// {
//     ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>(); 
//     CreateTessellation2D(tessellation);  
    
//     std::vector<double> x(tessellation->NbPoints()); 
//     for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)     x[i] = tessellation->Point(i).X(); 
    
//     return x; 
// }


// std::vector<double> EmbeddedIgaModeler::PrintNodesY()
// {
//     ANurbs::Pointer<ANurbs::CurveTessellation3D> tessellation = ANurbs::New<ANurbs::CurveTessellation3D>(); 
//     CreateTessellation2D(tessellation);  

//     std::vector<double> y(tessellation->NbPoints()); 
//     for (unsigned int i = 0; i < tessellation->NbPoints(); ++i)    y[i] = tessellation->Point(i).Y(); 
    
//     return y; 
// }

EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
      m_model_part(rModelPart)
{}
} // namespace Kratos.
