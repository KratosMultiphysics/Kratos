// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    JMCarbonell,
//                   Vahid Galavi
//

#pragma once

/* Project includes */
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
{
public:
    static double CalculateVonMisesStress(const Vector& StressVector);
    static double CalculateTrace(const Vector& StressVector);
    static double CalculateMeanStress(const Vector& StressVector);
    static double CalculateVonMisesStrain(const Vector& StrainVector);
    static Vector CalculateHenckyStrain(const Matrix& DeformationGradient, size_t VoigtSize);
    static Matrix CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient);
};

} // namespace Kratos