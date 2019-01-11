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
ANurbs::CurveTessellation3D EmbeddedIgaModeler::CreateTessellation()
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
                geometry->SetPole(node->GetId(),{node->X(), node->Y(), node->Z()});
            }
            KRATOS_WATCH()
             
             
             

            ANurbs::Curve3D curve(geometry); 
            const double flatness = 1e-2;
            ANurbs::CurveTessellation3D tessellation; 
            tessellation.Compute(curve, flatness); 
            
            
            KRATOS_WATCH(tessellation.Point(0).X())
            KRATOS_WATCH(tessellation.Point(0).Y())
            KRATOS_WATCH(tessellation.Point(0).Z())
            KRATOS_WATCH(tessellation.Point(1).X())
            KRATOS_WATCH(tessellation.Point(1).Y())
            KRATOS_WATCH(tessellation.Point(1).Z())
            
            return tessellation; 

        }
    }
}

void EmbeddedIgaModeler::CreateElements()
{
    // Perform Tessellation of Nurbs Curve
    auto tessellation = CreateTessellation(); 
    Model model; 
    ModelPart& skin_model_part = model.CreateModelPart("SkinModelPart"); 

    

    
    // // Create Nodes in Skin Modelpart from the Points generated in the tessellation of the curve    
    // for (int point_i = 0; point_i < tessellation.NbPoints(); ++point_i)    
    // {   
    //     skin_model_part.CreateNewNode(point_i, tessellation.Point(point_i).X(), tessellation.Point(point_i).Y(),tessellation.Point(point_i).Z()); 
    //     KRATOS_WATCH(skin_model_part.GetNode(point_i))
    // }
    // // skin_model_part property Container (Container is empty, but for the skin no properties are needed)
    // Properties::Pointer prop = skin_model_part.pGetProperties(0);
    
    // KRATOS_WATCH(tessellation.NbPoints())
    // int start_node_id = 0; 
    // int end_node_id = 1; 
    // // Create Elements in skin_model_part
    // for (int element_i = 0; element_i < tessellation.NbPoints() - 1; ++element_i)
    // {
    //     skin_model_part.CreateNewElement("Element2D2N", element_i, {{1, 2}}, prop); 
    //     start_node_id += 1; 
    //     end_node_id += 1; 

    //     KRATOS_WATCH(skin_model_part.GetElement(element_i))
    // }
    

}

EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
      m_model_part(rModelPart)
      

{
}
} // namespace Kratos.
