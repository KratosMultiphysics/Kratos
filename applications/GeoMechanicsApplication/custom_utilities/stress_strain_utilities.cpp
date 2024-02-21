// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "stress_strain_utilities.hpp"

namespace Kratos
{

Matrix StressStrainUtilities::CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient)
{
    return 0.5 * (prod(trans(rDeformationGradient), rDeformationGradient) -
                  IdentityMatrix(rDeformationGradient.size1()));
}

}
