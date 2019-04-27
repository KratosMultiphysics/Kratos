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
#include "embedded_iga_error_estimation.h"

namespace Kratos
{   
void EmbeddedIgaErrorEstimation::InsertGaussPointsExactSurface(
    const BrepFace& rFaceGeometry,
    const std::vector<Matrix>& rTriangulation_uv,
    std::vector<Matrix>& rGaussPoints_xyz)
{
    /**
     * This function inserts Gauss-Legendre Points into the triangulation of the surface in the parametric domain
     * and subsequently projects these points onto the NURBS surface.
    */

    const auto gp_canonical_tri = 
        Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

    
    std::vector<Matrix> gauss_points_uv(rTriangulation_uv.size(), ZeroMatrix(gp_canonical_tri.size(),2)); 

    for (unsigned int tri_i = 0; tri_i < rTriangulation_uv.size(); ++tri_i)
    {    
        for (unsigned int gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
        {
            gauss_points_uv[tri_i](gp_i,0) = rTriangulation_uv[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                             rTriangulation_uv[tri_i](1,0) * gp_canonical_tri[gp_i][0] + 
                                             rTriangulation_uv[tri_i](2,0) * gp_canonical_tri[gp_i][1];

            gauss_points_uv[tri_i](gp_i,1) = rTriangulation_uv[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                             rTriangulation_uv[tri_i](1,1) * gp_canonical_tri[gp_i][0] + 
                                             rTriangulation_uv[tri_i](2,1) * gp_canonical_tri[gp_i][1]; 
        }
    }

    EmbeddedIgaMapper::MapCartesianSpace(
        rFaceGeometry, gauss_points_uv, rGaussPoints_xyz); 
}

void EmbeddedIgaErrorEstimation::InsertGaussPointsApproxSurface(
    const BrepFace& rFaceGeometry,
    const std::vector<Matrix>& rTriangulation_uv,
    std::vector<Matrix>& rGaussPoints_xyz)
{
    /**
     * This function first maps the triangulation into the cartesian space and consequently inserts 
     * Gauss-Legendre points into the triangles approximating the exact Surface.
     * These points can be used to measure the distance between the approximated and exact surface
    */

    std::vector<Matrix> triangulation_xyz; 
    EmbeddedIgaMapper::MapCartesianSpace(
        rFaceGeometry, rTriangulation_uv, triangulation_xyz); 

    const auto gp_canonical_tri = 
        Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

    rGaussPoints_xyz.resize(triangulation_xyz.size(), ZeroMatrix(gp_canonical_tri.size(),3));  

    for (unsigned int tri_i = 0; tri_i < triangulation_xyz.size(); ++tri_i)
    {    
        for (unsigned int gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
        {
            rGaussPoints_xyz[tri_i](gp_i,0) = triangulation_xyz[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                              triangulation_xyz[tri_i](1,0) * gp_canonical_tri[gp_i][0] + 
                                              triangulation_xyz[tri_i](2,0) * gp_canonical_tri[gp_i][1];

            rGaussPoints_xyz[tri_i](gp_i,1) = triangulation_xyz[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                              triangulation_xyz[tri_i](1,1) * gp_canonical_tri[gp_i][0] + 
                                              triangulation_xyz[tri_i](2,1) * gp_canonical_tri[gp_i][1]; 
        
            rGaussPoints_xyz[tri_i](gp_i,2) = triangulation_xyz[tri_i](0,2) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                              triangulation_xyz[tri_i](1,2) * gp_canonical_tri[gp_i][0] + 
                                              triangulation_xyz[tri_i](2,2) * gp_canonical_tri[gp_i][1]; 
        }
    }
}

void EmbeddedIgaErrorEstimation::GetError(
    const std::vector<Matrix>& rGaussPointsExact, 
    const std::vector<Matrix>& rGaussPointsApprox, 
    Vector& rError)
{   
    const auto number_tri = rGaussPointsExact.size();
    const auto number_points = rGaussPointsExact[0].size1();
    rError = ZeroVector(number_tri); 

    double ele_error; 

    for (unsigned int tri_i = 0; tri_i < number_tri; ++tri_i)
    {
        ele_error = 0; 
        for (unsigned int point_i = 0; point_i < number_points; ++point_i)
        {
            ele_error += sqrt(pow((rGaussPointsExact[tri_i](point_i, 0) - rGaussPointsApprox[tri_i](point_i, 0)),2) + 
                              pow((rGaussPointsExact[tri_i](point_i, 1) - rGaussPointsApprox[tri_i](point_i, 1)),2) + 
                              pow((rGaussPointsExact[tri_i](point_i, 2) - rGaussPointsApprox[tri_i](point_i, 2)),2)); 
        }
        rError[tri_i] = ele_error; 
    }
}

    EmbeddedIgaErrorEstimation::EmbeddedIgaErrorEstimation()
    {}

} // namespace Kratos.
