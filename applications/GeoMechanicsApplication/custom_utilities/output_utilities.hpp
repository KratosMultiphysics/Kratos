// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/kratos_export_api.h"
#include "includes/properties.h"

#include <algorithm>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) OutputUtilities
{
public:
    template <typename StressVectorContainer, typename OutputIter>
    static void CalculateShearCapacityValues(const StressVectorContainer& rStressVectors,
                                             OutputIter                   Destination,
                                             const Properties&            rProperties)
    {
        const auto c   = ConstitutiveLawUtilities::GetCohesion(rProperties);
        const auto phi = ConstitutiveLawUtilities::GetFrictionAngleInDegrees(rProperties);
        auto       calculate_shear_capacity = [c, phi](const auto& rStressVector) {
            return StressStrainUtilities::CalculateMohrCoulombShearCapacity(rStressVector, c, phi);
        };
        std::transform(std::begin(rStressVectors), std::end(rStressVectors), Destination, calculate_shear_capacity);
    }
};

} // namespace Kratos