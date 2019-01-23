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
    void EmbeddedIgaErrorEstimation::InsertGaussPoints(std::vector<array_1d<double, 2> >& gp_pos)
    {
        auto gp_coords = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        gp_pos.resize(mTriangles.size() * gp_coords.size()); 
        
        unsigned int id = 0; 
        for (unsigned int tri_i = 0; tri_i < mTriangles.size(); ++tri_i)
        {    
            for (unsigned int gp_i = 0; gp_i < gp_coords.size(); ++gp_i)
            {
                gp_pos[id][0] = mTriangles[tri_i](0,0) * (1 - gp_coords[gp_i][0] - gp_coords[gp_i][1]) + mTriangles[tri_i](1,0) * gp_coords[gp_i][0] + mTriangles[tri_i](2,0) * gp_coords[gp_i][1]; 
                gp_pos[id][1] = mTriangles[tri_i](0,1) * (1 - gp_coords[gp_i][0] - gp_coords[gp_i][1]) + mTriangles[tri_i](1,1) * gp_coords[gp_i][0] + mTriangles[tri_i](2,1) * gp_coords[gp_i][1]; 
                id += 1; 
            }
        }
    }

    EmbeddedIgaErrorEstimation::EmbeddedIgaErrorEstimation(std::vector<Matrix> rTriangles)
        : mTriangles(rTriangles)
    {
    }

} // namespace Kratos.
