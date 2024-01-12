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

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ExpectCouplingTermsWhenNotCompatibleWithPlaxis, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy");

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
    StructuredMeshGeneratorProcess(domain_geometry, r_model_part, mesher_parameters).Execute();

    auto p_dummy_law              = std::make_shared<StubLinearElasticLaw>();
    auto& r_model_part_properties = r_model_part.GetProperties(0);
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);

    Parameters k0_settings;

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    const auto& r_elements = r_model_part.Elements();
    KRATOS_EXPECT_TRUE(std::all_of(r_elements.begin(), r_elements.end(), [](const auto& rElement) {
        auto p_constitutive_law =
            dynamic_cast<GeoLinearElasticLaw*>(rElement.GetProperties().GetValue(CONSTITUTIVE_LAW).get());
        return p_constitutive_law->GetCouplingBehavior() == GeoLinearElasticLaw::Coupling::Yes;
    }))
}
} // namespace Kratos::Testing
