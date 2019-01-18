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
                    
                    rPolygon.resize(point_id + tessellation->NbPoints() - 1); 

                    KRATOS_WATCH(tessellation->NbPoints())
                    
                    for (unsigned int i = 0; i < tessellation->NbPoints() -1; ++i)
                    {
                        rPolygon[point_id][0] = tessellation->Point(i).X(); 
                        rPolygon[point_id][1] = tessellation->Point(i).Y(); 
                        // rPolygon[point_id][2] = tessellation->Point(i).Z();
                        point_id += 1; 
                    }
                }
            }
        }
    }

}


bool EmbeddedIgaModeler::CreateElements3D(const std::vector<array_1d<double,3>>& rPolygon, std::vector<Matrix>& rTriangles)
{
    array_1d<double,3> p1,p2,p3,p4; 
    int bestvertex;
    double weight, minweight, d1, d2;
    Diagonal diagonal, newdiagonal;
    std::list<Diagonal> diagonals;
    bool ret = true;
    const unsigned int number_points = rPolygon.size();
    matrix<DPState> dpstates(number_points, number_points);

    for (unsigned int i = 0; i < number_points; ++i)
    {   
        p1 = rPolygon[i];  
        for(unsigned int j = i + 1; j < number_points; ++j)
        {
            dpstates(j,i).visible = true;
            dpstates(j,i).weight = 0;
            dpstates(j,i).bestvertex = -1;

            if (j != i+1)
            {
                p2 = rPolygon[j]; 
                if (i == 0)     p3 = rPolygon[number_points-1]; 
                else            p3 = rPolygon[i - 1]; 
                
                if (i == number_points - 1)     p4 = rPolygon[0]; 
                else                            p4 = rPolygon[i + 1]; 
                
                if (!InCone(p3, p1, p4, p2))
                {           
                    dpstates(j,i).visible = false; 
                    continue; 
                }

                if (j == 0)     p3 = rPolygon[number_points - 1]; 
                else            p3 = rPolygon[j - 1]; 

                if (j == (number_points - 1))       p4 = rPolygon[0];
                else                                p4 = rPolygon[j + 1];
                
                if (!InCone(p3, p2, p4, p1)) 
                {
                    dpstates(j, i).visible = false;
                    continue;
                }

                for (unsigned int k = 0; k < number_points; ++k) 
                {
                    p3 = rPolygon[k];
                    if (k == number_points - 1)   p4 = rPolygon[0];
                    else                            p4 = rPolygon[k + 1];

                    if (Intersects(p1, p2, p3, p4)) 
                    {
                        dpstates(j, i).visible = false;
                        break;
                    }
                }
            }
        }      
    }

    dpstates(number_points - 1, 0).visible = true;
    dpstates(number_points - 1, 0).weight = 0;
    dpstates(number_points - 1, 0).bestvertex = -1;

    for (unsigned int gap = 2; gap < number_points; ++gap) 
    {
        for (unsigned int i = 0; i < number_points - gap; ++i) 
        {
            int j = i + gap;
            
            if (!dpstates(j, i).visible)    continue;

            int bestvertex = -1;
        
            for (unsigned int k = (i + 1); k<j; k++) 
            {
                if (!dpstates(k, i).visible)    continue;
                if (!dpstates(j, k).visible)    continue;

                if (k <= i + 1)     d1 = 0;
                else                d1 = Distance(rPolygon[i], rPolygon[k]);

                if (j <= k + 1)     d2 = 0;
                else d2 = Distance(rPolygon[k], rPolygon[j]);

                weight = dpstates(k,i).weight + dpstates(j,k).weight + d1 + d2;


                if (bestvertex == -1 || weight < minweight) 
                {
                    bestvertex = k;
                    minweight = weight;
                }
            }
            if (bestvertex == -1)       return false; 
            
            dpstates(j,i).bestvertex = bestvertex;
            dpstates(j,i).weight = minweight;
        }
    }

    newdiagonal.index1 = 0;
    newdiagonal.index2 = number_points - 1;
    diagonals.push_back(newdiagonal);

    while (!diagonals.empty()) 
    {
        diagonal = *(diagonals.begin());
        diagonals.pop_front();
        bestvertex = dpstates(diagonal.index2, diagonal.index1).bestvertex;
        if (bestvertex == -1) 
        {
        ret = false;
        break;
        }

        Matrix triangle(3, 2);
        triangle(0, 0) = rPolygon[diagonal.index1][0];
        triangle(0, 1) = rPolygon[diagonal.index1][1];
        triangle(1, 0) = rPolygon[bestvertex][0];
        triangle(1, 1) = rPolygon[bestvertex][1];
        triangle(2, 0) = rPolygon[diagonal.index2][0];
        triangle(2, 1) = rPolygon[diagonal.index2][1];
 
        if (abs(GetAreaOfTriangle(triangle))>1e-9)
        {
            rTriangles.push_back(triangle);
        }
        else
        {
            std::cout << "triangle with zero area" << GetAreaOfTriangle(triangle) << std::endl;
            KRATOS_WATCH(triangle)
        }
        if (bestvertex > (diagonal.index1 + 1)) 
        {
            newdiagonal.index1 = diagonal.index1;
            newdiagonal.index2 = bestvertex;
            diagonals.push_back(newdiagonal);
        }

        if (diagonal.index2 > (bestvertex + 1)) 
        {
            newdiagonal.index1 = bestvertex;
            newdiagonal.index2 = diagonal.index2;
            diagonals.push_back(newdiagonal);
        }
    }
    
    return true;
}


std::vector<Matrix> EmbeddedIgaModeler::Triangulate()
{
    std::vector<array_1d<double,3>> polygon; 
    CreateTessellationParameterCurve(polygon); 
    std::vector<Matrix> triangles;
    bool success = CreateElements3D(polygon, triangles);

    KRATOS_WATCH(triangles)
    KRATOS_WATCH(triangles.size())
    
    return triangles;
}

bool EmbeddedIgaModeler::IsConvex(
    const array_1d<double, 3>& p1, const array_1d<double, 3>& p2, 
    const array_1d<double, 3>& p3)
{   
    double tmp = (p3[1] - p1[1])*(p2[0] - p1[0]) - (p3[0] - p1[0])*(p2[1] - p1[1]);

    if (tmp > 0)       return true;
    else                return false;
}

bool EmbeddedIgaModeler::InCone(
    array_1d<double, 3> &p1, array_1d<double, 3> &p2,
    array_1d<double, 3> &p3, array_1d<double, 3> &p) 
{
    if (IsConvex(p1, p2, p3)) {
      if (!IsConvex(p1, p2, p)) return false;
      if (!IsConvex(p2, p3, p)) return false;
      return true;
    }
    else {
      if (IsConvex(p1, p2, p)) return true;
      if (IsConvex(p2, p3, p)) return true;
      return false;
    }
}

bool EmbeddedIgaModeler::Intersects(
    array_1d<double, 3>& p11, array_1d<double, 3>& p12,
    array_1d<double, 3>& p21, array_1d<double, 3>& p22)
{
    if ((p11[0] == p21[0]) && (p11[1] == p21[1])) return false;
    if ((p11[0] == p22[0]) && (p11[1] == p22[1])) return false;
    if ((p12[0] == p21[0]) && (p12[1] == p21[1])) return false;
    if ((p12[0] == p22[0]) && (p12[1] == p22[1])) return false;

    array_1d<double, 2> v1ort, v2ort, v;
    double dot11, dot12, dot21, dot22;

    v1ort[0] = p12[1] - p11[1];
    v1ort[1] = p11[0] - p12[0];

    v2ort[0] = p22[1] - p21[1];
    v2ort[1] = p21[0] - p22[0];

    v[0] = p21[0] - p11[0];
    v[1] = p21[1] - p11[1];
    dot21 = v[0] * v1ort[0] + v[1] * v1ort[1];
    v[0] = p22[0] - p11[0];
    v[1] = p22[1] - p11[1];
    dot22 = v[0] * v1ort[0] + v[1] * v1ort[1];

    v[0] = p11[0] - p21[0];
    v[1] = p11[1] - p21[1];
    dot11 = v[0] * v2ort[0] + v[1] * v2ort[1];
    v[0] = p12[0] - p21[0];
    v[1] = p12[1] - p21[1];
    dot12 = v[0] * v2ort[0] + v[1] * v2ort[1];

    if (dot11*dot12>0) return false;
    if (dot21*dot22>0) return false;

    return true;
}

double EmbeddedIgaModeler::Distance(array_1d<double, 3> p1, array_1d<double, 3> p2)
{
    return sqrt(p1[0] * p2[0] + p1[1] * p2[1]);
}

double EmbeddedIgaModeler::GetAreaOfTriangle(const Matrix& triangle)
  {
    double area = abs((triangle(0, 0)*(triangle(1, 1) - triangle(2, 1))
                + triangle(1, 0)*(triangle(2, 1) - triangle(0, 1))
                + triangle(2, 0)*(triangle(0, 1) - triangle(1, 1))) / 2);

    return area; 
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
    std::vector<Matrix> triangles = Triangulate();
    std::vector<double> x(triangles.size()*3);
    int index = 0; 
    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            x[index] = triangles[i](j,0);
            index += 1; 
        }
    }
    return x; 
}

std::vector<double> EmbeddedIgaModeler::PrintNodesY3D()
{
    std::vector<Matrix> triangles = Triangulate();
    std::vector<double> y(triangles.size()*3);
    int index = 0; 
    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            y[index] = triangles[i](j,1);
            index += 1; 
        }
    }
    return y; 
}

std::vector<double> EmbeddedIgaModeler::TessellationX()
{
    std::vector<array_1d<double,3>> polygon; 
    CreateTessellationParameterCurve(polygon); 
    std::vector<Matrix> triangles;
    bool success = CreateElements3D(polygon, triangles);
    std::vector<double> x(polygon.size()); 
    KRATOS_WATCH(polygon.size())

    for (int i = 0; i < polygon.size(); ++i)
    {
        x[i] = polygon[i][0]; 
    }

    return x; 
}

std::vector<double> EmbeddedIgaModeler::TessellationY()
{
    std::vector<array_1d<double,3>> polygon; 
    CreateTessellationParameterCurve(polygon); 
    std::vector<Matrix> triangles;
    bool success = CreateElements3D(polygon, triangles);
    std::vector<double> y(polygon.size()); 
    KRATOS_WATCH(polygon.size())

    for (int i = 0; i < polygon.size(); ++i)
    {
        y[i] = polygon[i][0]; 
    }

    return y; 
}




EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
      m_model_part(rModelPart)
{
}
} // namespace Kratos.
