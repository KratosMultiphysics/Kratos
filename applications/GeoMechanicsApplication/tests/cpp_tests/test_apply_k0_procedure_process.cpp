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
#include "containers/model.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "stub_linear_elastic_law.h"
#include "testing/testing.h"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>

using namespace Kratos;

namespace
{

ModelPart& PrepareTestModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy");

    // Set up the test model part mesh
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, result, mesher_parameters).Execute();

    auto p_dummy_law              = std::make_shared<Testing::StubLinearElasticLaw>();
    auto& r_model_part_properties = result.GetProperties(0);
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);

    return result;
}

bool ElementConsidersDiagonalEntriesOnlyAndNoShear(const Element& rElement)
{
    auto p_constitutive_law =
        dynamic_cast<const GeoLinearElasticLaw*>(rElement.GetProperties().GetValue(CONSTITUTIVE_LAW).get());
    return p_constitutive_law && p_constitutive_law->GetConsiderDiagonalEntriesOnlyAndNoShear();
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AllElementsConsiderDiagonalEntriesOnlyAndNoShearWhenUseStandardProcedureFlagIsNotDefined,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings; // 'use_standard_procedure' is not defined, assume it to be true

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::all_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(AllElementsConsiderDiagonalEntriesOnlyAndNoShearWhenUsingStandardProcedure,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::all_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(NoneOfElementsConsiderDiagonalEntriesOnlyAndNoShearWhenNotUsingStandardProcedure,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": false})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::none_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(UseStandardProcedureFlagIsInEffectDuringProcessExecutionOnly, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize(); // start considering diagonal entries only and no shear
    process.ExecuteFinalize(); // stop considering diagonal entries only and no shear

    KRATOS_EXPECT_TRUE(boost::algorithm::none_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

} // namespace Kratos::Testing
