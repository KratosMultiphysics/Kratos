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
#include "embedded_iga_triangulation.h"

namespace Kratos
{   

    bool EmbeddedIgaTriangulation::CreateTriangles(const std::vector<array_1d<double,3>>& rPolygon, std::vector<Matrix>& rTriangles)
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

                if (j != i + 1)
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

            // Matrix with the triangle coordinates
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

    bool EmbeddedIgaTriangulation::IsConvex(
        const array_1d<double, 3>& p1, const array_1d<double, 3>& p2, 
        const array_1d<double, 3>& p3)
    {   
        double tmp = (p3[1] - p1[1])*(p2[0] - p1[0]) - (p3[0] - p1[0])*(p2[1] - p1[1]);

        if (tmp > 0)    return true;
        else            return false;
    }

    bool EmbeddedIgaTriangulation::InCone(
        array_1d<double, 3>& p1, array_1d<double, 3>& p2,
        array_1d<double, 3>& p3, array_1d<double, 3>& p) 
    {
        if (IsConvex(p1, p2, p3)) 
        {
            if (!IsConvex(p1, p2, p)) return false;
            if (!IsConvex(p2, p3, p)) return false;

            return true;
        }
        else 
        {
            if (IsConvex(p1, p2, p)) return true;
            if (IsConvex(p2, p3, p)) return true;
            
            return false;
        }
    }

    bool EmbeddedIgaTriangulation::Intersects(
        array_1d<double, 3>& p11, array_1d<double, 3>& p12,
        array_1d<double, 3>& p21, array_1d<double, 3>& p22)
    {
        if (p11[0] == p21[0] && p11[1] == p21[1]) return false;
        if (p11[0] == p22[0] && p11[1] == p22[1]) return false;
        if (p12[0] == p21[0] && p12[1] == p21[1]) return false;
        if (p12[0] == p22[0] && p12[1] == p22[1]) return false;

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

        if (dot11 * dot12 > 0) return false;
        if (dot21 * dot22 > 0) return false;

        return true;
    }

    double EmbeddedIgaTriangulation::Distance(array_1d<double, 3> p1, array_1d<double, 3> p2)
    {
        return sqrt(p1[0] * p2[0] + p1[1] * p2[1]);
    }

    double EmbeddedIgaTriangulation::GetAreaOfTriangle(const Matrix& triangle)
    {
        double area = abs((triangle(0, 0)*(triangle(1, 1) - triangle(2, 1))
                    + triangle(1, 0)*(triangle(2, 1) - triangle(0, 1))
                    + triangle(2, 0)*(triangle(0, 1) - triangle(1, 1))) / 2);

        return area; 
    }

} // namespace Kratos.
