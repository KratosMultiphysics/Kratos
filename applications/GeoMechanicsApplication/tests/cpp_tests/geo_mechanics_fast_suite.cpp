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
//                   Richard Faasse
//

#include "geo_mechanics_fast_suite.h"
#include "geo_mechanics_application.h"
#include "linear_solvers_application.h"

namespace Kratos::Testing
{

KratosGeoMechanicsFastSuite::KratosGeoMechanicsFastSuite() : KratosCoreFastSuite()
{
    mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpGeoApp);
    mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
    this->ImportApplicationIntoKernel(mpLinearSolversApp);
}

KratosGeoMechanicsFastSuiteWithoutKernel::KratosGeoMechanicsFastSuiteWithoutKernel()
    : KratosCoreFastSuiteWithoutKernel()
{
    // clang-format off
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_ACCELERATION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TOTAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INCREMENTAL_DISPLACEMENT)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TOTAL_ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INCREMENTAL_ROTATION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FLUID_FLUX_VECTOR)
    KRATOS_REGISTER_VARIABLE(HYDRAULIC_HEAD)
    KRATOS_REGISTER_VARIABLE(WATER_PRESSURE)

    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(TIME_STEPS)

    KRATOS_REGISTER_VARIABLE(GEO_COULOMB_HARDENING_TYPE)
    KRATOS_REGISTER_VARIABLE(GEO_FRICTION_ANGLE)
    KRATOS_REGISTER_VARIABLE(GEO_COHESION)
    KRATOS_REGISTER_VARIABLE(GEO_DILATANCY_ANGLE)
    KRATOS_REGISTER_VARIABLE(GEO_TENSILE_STRENGTH)
    KRATOS_REGISTER_VARIABLE(GEO_COULOMB_HARDENING_MAX_ITERATIONS)
    KRATOS_REGISTER_VARIABLE(GEO_COULOMB_HARDENING_CONVERGENCE_TOLERANCE)
    // clang-format on
}

KratosGeoMechanicsIntegrationSuite::KratosGeoMechanicsIntegrationSuite()
    : KratosCoreFastSuite()
{
    mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpGeoApp);
    mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
    this->ImportApplicationIntoKernel(mpLinearSolversApp);
}

} // namespace Kratos::Testing
