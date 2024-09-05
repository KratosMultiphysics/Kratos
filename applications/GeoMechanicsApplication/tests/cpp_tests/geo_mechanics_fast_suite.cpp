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

KratosGeoMechanicsFastSuiteWithoutKernel::KratosGeoMechanicsFastSuiteWithoutKernel() : KratosCoreFastSuiteWithoutKernel()
{
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_ACCELERATION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FLUID_FLUX_VECTOR)
    KRATOS_REGISTER_VARIABLE(HYDRAULIC_HEAD)

    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_TENSOR)
}

KratosGeoMechanicsIntegrationSuite::KratosGeoMechanicsIntegrationSuite() : KratosCoreFastSuite()
{
    mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpGeoApp);
    mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
    this->ImportApplicationIntoKernel(mpLinearSolversApp);
}

} // namespace Kratos::Testing
