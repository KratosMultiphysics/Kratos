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


    std::vector<std::vector<double>> EmbeddedIgaModeler::TestTriangle()
    // void EmbeddedIgaModeler::TestTriangle()
    {
        std::vector<array_1d<double,3>> polygon;
        CreateTessellationParameterCurve(polygon);

        // initializing the i/o containers
        struct triangulateio in_data; 
        struct triangulateio out_data; 
        struct triangulateio vor_out_data;

        InitTriangulationDataStructure(in_data); 
        InitTriangulationDataStructure(out_data); 
        InitTriangulationDataStructure(vor_out_data); 

        // Initialize the pointlist (1d list) with the number of points and the position 
        in_data.numberofpoints = polygon.size(); 
        in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));

        unsigned int point_idx = 0;
        for (unsigned int i = 0; i < in_data.numberofpoints; ++i)
        {
            for (unsigned int j = 0; j < 2; ++j)    in_data.pointlist[point_idx++] = polygon[i][j];
        }
        
        // Initilize the segment list with the number of boundary edges and the start and end node id
        // For closed polygons the number of segments is equal to the number of points
        in_data.numberofsegments = in_data.numberofpoints; 
        in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        
        unsigned int vertex_id = 1;
        unsigned int seg_idx = 1;  
        in_data.segmentlist[0] = 0; 
        for (unsigned int seg_idx = 1; seg_idx < in_data.numberofsegments * 2 - 1;)
        {
            for (unsigned int j = 0; j < 2; ++j)    in_data.segmentlist[seg_idx++] =  vertex_id;
            vertex_id += 1; 
        }
        in_data.segmentlist[in_data.numberofsegments * 2 - 1] = 0; 

        char trigenOptsVerbose[] = "pz"; 
        char* trigenOpts = trigenOptsVerbose; 

        triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

        std::vector<std::vector<double>> tri_coords (out_data.numberoftriangles * 3 , std::vector<double>(2,0)); 
        unsigned int id = 0; 
        for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                tri_coords[id + j][0] = out_data.pointlist[out_data.trianglelist[id + j] * 2];
                tri_coords[id + j][1] = out_data.pointlist[out_data.trianglelist[id + j] * 2 + 1]; 
            }
            id += 3; 
        }
        return tri_coords; 
    }

    




    void EmbeddedIgaModeler::InitTriangulationDataStructure(triangulateio& tr)
    {
        tr.pointlist                  = (REAL*) NULL;
        tr.pointattributelist         = (REAL*) NULL;
        tr.pointmarkerlist            = (int*) NULL;
        tr.numberofpoints             = 0;
        tr.numberofpointattributes    = 0;
        tr.trianglelist               = (int*) NULL;
        tr.triangleattributelist      = (REAL*) NULL;
        tr.trianglearealist           = (REAL*) NULL;
        tr.neighborlist               = (int*) NULL;
        tr.numberoftriangles          = 0;
        tr.numberofcorners            = 3;
        tr.numberoftriangleattributes = 0;
        tr.segmentlist                = (int*) NULL;
        tr.segmentmarkerlist          = (int*) NULL;
        tr.numberofsegments           = 0;
        tr.holelist                   = (REAL*) NULL;
        tr.numberofholes              = 0;
        tr.regionlist                 = (REAL*) NULL;
        tr.numberofregions            = 0;
        tr.edgelist                   = (int*) NULL;
        tr.edgemarkerlist             = (int*) NULL;
        tr.normlist                   = (REAL*) NULL;
        tr.numberofedges              = 0;
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