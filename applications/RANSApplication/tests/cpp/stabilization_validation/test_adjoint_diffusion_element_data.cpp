//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/properties.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/adjoint_test_utilities.h"
#include "custom_utilities/fluid_test_utilities.h"
#include "rans_application_variables.h"

// Element data containers
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data_derivatives.h"

namespace Kratos
{
namespace Testing
{

namespace
{
void StabilizationValidationAddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
}

void StabilizationValidationSetVariables(ModelPart& rModelPart)
{
    using namespace RansApplicationTestUtilities;

    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY_POTENTIAL, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
}

void StabilizationValidationSetProperties(Properties& rProperties)
{
    rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("Newtonian2DLaw").Clone());
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(StabilizationValidationDiffusionPhiDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        StabilizationValidationElementData::DiffusionElementData, StabilizationValidationElementData::DiffusionElementDataDerivatives::PhiDerivative>(
        model, StabilizationValidationAddVariables, StabilizationValidationSetVariables,
        StabilizationValidationSetProperties, [](ModelPart&){}, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(StabilizationValidationDiffusionShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        StabilizationValidationElementData::DiffusionElementData, StabilizationValidationElementData::DiffusionElementDataDerivatives::ShapeDerivative>(
        model, StabilizationValidationAddVariables, StabilizationValidationSetVariables,
        StabilizationValidationSetProperties, [](ModelPart&){}, 2, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos