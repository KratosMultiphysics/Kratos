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
    void EmbeddedIgaErrorEstimation::InsertGaussPoints(std::vector<Matrix>& rGP_uv)
    {
        const auto gp_canonical_tri = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        rGP_uv.resize(mTriangles.size(), ZeroMatrix(gp_canonical_tri.size(),2));  
        
        for (unsigned int tri_i = 0; tri_i < mTriangles.size(); ++tri_i)
        {    
            for (unsigned int gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
            {
                rGP_uv[tri_i](gp_i,0) =  mTriangles[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                         mTriangles[tri_i](1,0) * gp_canonical_tri[gp_i][0] + 
                                         mTriangles[tri_i](2,0) * gp_canonical_tri[gp_i][1];

                rGP_uv[tri_i](gp_i,1) =  mTriangles[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) + 
                                         mTriangles[tri_i](1,1) * gp_canonical_tri[gp_i][0] + 
                                         mTriangles[tri_i](2,1) * gp_canonical_tri[gp_i][1]; 
            }
        }
    }

    EmbeddedIgaErrorEstimation::EmbeddedIgaErrorEstimation(std::vector<Matrix> rTriangles)
        : mTriangles(rTriangles)
    {
    }

} // namespace Kratos.
