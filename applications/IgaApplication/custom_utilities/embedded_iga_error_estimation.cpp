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
    void EmbeddedIgaErrorEstimation::InsertGaussPointsSurface(std::vector<Matrix>& rGP_surface)
    {
        /**
         * This function inserts Gauss-Legendre Points on the exact surface described by NURBS.
         * This is achieved by inserting the points on the triangulated parameter space and
         * projecting theses points onto the surface. 
         * 
         * => TODO: Projection of the points onto the surface
        */

        const auto gp_canonical_tri = 
            Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        rGP_surface.resize(mTriangles.size(), ZeroMatrix(gp_canonical_tri.size(),2)); 

        KRATOS_WATCH(gp_canonical_tri.size())
        
        for (unsigned int tri_i = 0; tri_i < mTriangles.size(); ++tri_i)
        {    
            for (unsigned int gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
            {
                rGP_surface[tri_i](gp_i,0) = mTriangles[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                             mTriangles[tri_i](1,0) * gp_canonical_tri[gp_i][0] + 
                                             mTriangles[tri_i](2,0) * gp_canonical_tri[gp_i][1];

                rGP_surface[tri_i](gp_i,1) = mTriangles[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                             mTriangles[tri_i](1,1) * gp_canonical_tri[gp_i][0] + 
                                             mTriangles[tri_i](2,1) * gp_canonical_tri[gp_i][1]; 
            }
        }
    }

    void EmbeddedIgaErrorEstimation::InsertGaussPointsTriangulation(std::vector<Matrix>& rGP_tri)
    {
        /**
         * This function inserts Gauss-Legendre points into the triangles approximating the exact Surface.
         * These points can be used to measure the distance between the approximated and exact surface
        */

        const auto gp_canonical_tri = 
            Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        rGP_tri.resize(mTriangles.size(), ZeroMatrix(gp_canonical_tri.size(),3));  

        for (unsigned int tri_i = 0; tri_i < mTriangles.size(); ++tri_i)
        {    
            for (unsigned int gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
            {
                rGP_tri[tri_i](gp_i,0) = mTriangles[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                         mTriangles[tri_i](1,0) * gp_canonical_tri[gp_i][0] + 
                                         mTriangles[tri_i](2,0) * gp_canonical_tri[gp_i][1];

                rGP_tri[tri_i](gp_i,1) = mTriangles[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                         mTriangles[tri_i](1,1) * gp_canonical_tri[gp_i][0] + 
                                         mTriangles[tri_i](2,1) * gp_canonical_tri[gp_i][1]; 
            
                rGP_tri[tri_i](gp_i,2) = mTriangles[tri_i](0,2) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                         mTriangles[tri_i](1,2) * gp_canonical_tri[gp_i][0] + 
                                         mTriangles[tri_i](2,2) * gp_canonical_tri[gp_i][1]; 
            }
        }
    }




    EmbeddedIgaErrorEstimation::EmbeddedIgaErrorEstimation(std::vector<Matrix> rTriangles)
        : mTriangles(rTriangles)
    {
    }

} // namespace Kratos.
