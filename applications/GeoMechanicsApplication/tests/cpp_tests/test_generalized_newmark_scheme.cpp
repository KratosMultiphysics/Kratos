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

#include "custom_strategies/schemes/generalized_newmark_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing
{

// We need this class to test all non-abstract functions of the
// GeneralizedNewmarkScheme class. We cannot use the
// GeneralizedNewmarkScheme class directly, because it is abstract.
class ConcreteGeneralizedNewmarkScheme : public GeneralizedNewmarkScheme<SparseSpaceType, LocalSpaceType>
{
public:
    using GeneralizedNewmarkScheme::GeneralizedNewmarkScheme;

protected:
    void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        // Intentionally left empty
    }
};

KRATOS_TEST_CASE_IN_SUITE(ForInvalidTheta_CheckNewmarkScheme_Throws, KratosGeoMechanicsFastSuite)
{
    constexpr double invalid_theta = -2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConcreteGeneralizedNewmarkScheme scheme({}, invalid_theta),
                                      "Theta must be larger than zero, but got -2")
}

} // namespace Kratos::Testing
