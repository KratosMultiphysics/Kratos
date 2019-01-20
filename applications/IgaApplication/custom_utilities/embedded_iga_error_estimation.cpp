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
    void EmbeddedIgaErrorEstimation::PrintTriangleGaussPoints()
    {

        // xiTri = x1*(1-triGP(iGP,2)-triGP(iGP,3)) + x2 * triGP(iGP,2) + x3 * triGP(iGP,3); 
        // etaTri = y1*(1-triGP(iGP,2)-triGP(iGP,3)) + y2 * triGP(iGP,2) + y3 * triGP(iGP,3);

        



        // auto test = Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
        
        
    }

    EmbeddedIgaErrorEstimation::EmbeddedIgaErrorEstimation(std::vector<Matrix> rTriangles)
        : mTriangles(rTriangles)
    {
    }

} // namespace Kratos.
