//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// Project includes
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "modeler/gap_surrogate_boundary_modeler.h"
#include "utilities/gap_construction_utility.h"

namespace Kratos {

namespace {

Parameters MakeBaseModelerParameters(const Parameters& rParameters)
{
    Parameters base_parameters = rParameters.Clone();
    const std::array<std::string, 5> gap_only_entries{{
        "gap_element_sub_model_part_name",
        "gap_interface_sub_model_part_name",
        "gap_element_name",
        "gap_interface_condition_name",
        "gap_integration_order"}};
    for (const auto& r_entry : gap_only_entries) {
        base_parameters.RemoveValue(r_entry);
    }

    // Endpoint projections are the four corners of the ruled gap patch. This
    // internal setting is forced by the derived modeler, irrespective of a
    // value supplied by the caller.
    if (base_parameters.Has("project_surrogate_face_end_points")) {
        base_parameters["project_surrogate_face_end_points"].SetBool(true);
    } else {
        base_parameters.AddBool("project_surrogate_face_end_points", true);
    }
    return base_parameters;
}

} // unnamed namespace

GapSurrogateBoundaryModeler::GapSurrogateBoundaryModeler(
    Model& rModel,
    Parameters ModelerParameters)
    : SurrogateBoundaryModeler(
          rModel, MakeBaseModelerParameters(ModelerParameters)),
      mpGapModel(&rModel),
      mGapParameters(ModelerParameters)
{
    mGapParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

const Parameters GapSurrogateBoundaryModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "echo_level" : 0,
        "domain_size" : 2,
        "input_model_part_name" : "SkinModelPart",
        "output_model_part_name" : "SurrogateModelPart",
        "volume_sub_model_part_name" : "FluidDomain",
        "surrogate_sub_model_part_name" : "SurrogateBoundary",
        "lower_point" : [0.0, 0.0, 0.0],
        "upper_point" : [1.0, 1.0, 0.0],
        "number_of_elements" : [10, 10, 1],
        "physical_domain" : "outside",
        "lambda" : 1.0,
        "classification_samples_per_direction" : 3,
        "element_name" : "Element2D4N",
        "condition_name" : "LineCondition2D2N",
        "properties_id" : 1,
        "number_of_boundary_gauss_points" : 2,
        "project_surrogate_face_end_points" : true,
        "is_moving" : false,
        "background_sub_model_part_name" : "BackgroundGrid",
        "gap_element_sub_model_part_name" : "GapElements",
        "gap_interface_sub_model_part_name" : "GapInterfaces",
        "gap_element_name" : "",
        "gap_interface_condition_name" : "",
        "gap_integration_order" : 2
    })");
}

void GapSurrogateBoundaryModeler::SetupModelPart()
{
    KRATOS_ERROR_IF(mpGapModel == nullptr)
        << "GapSurrogateBoundaryModeler has no Model container." << std::endl;
    KRATOS_ERROR_IF(mGapParameters["domain_size"].GetInt() != 2)
        << "GapSurrogateBoundaryModeler currently supports domain_size=2 only."
        << std::endl;
    KRATOS_ERROR_IF(mGapParameters["gap_integration_order"].GetInt() != 2)
        << "GapSurrogateBoundaryModeler currently supports a 2x2 gap "
        << "quadrature only." << std::endl;
    KRATOS_ERROR_IF(mGapParameters["gap_element_name"].GetString().empty())
        << "GapSurrogateBoundaryModeler requires 'gap_element_name'."
        << std::endl;

    SurrogateBoundaryModeler::SetupModelPart();
    SetupGapModelPart2D();
}

void GapSurrogateBoundaryModeler::PrepareGeometryModel()
{
    // The standard Kratos lifecycle calls this stage once before
    // SetupModelPart. Match the base-modeler no-op on that initial call.
    if (!IsMovingModelInitialized()) {
        return;
    }
    KRATOS_ERROR_IF_NOT(mGapParameters["is_moving"].GetBool())
        << "GapSurrogateBoundaryModeler::PrepareGeometryModel requires "
        << "'is_moving' to be true." << std::endl;

    // Gap entities hold pointers to the old surrogate conditions and active
    // parent elements. Destroy this transient layer before the base modeler
    // replaces the active view and its boundary.
    ClearGapModelPart2D();
    SurrogateBoundaryModeler::PrepareGeometryModel();
    SetupGapModelPart2D();
}

void GapSurrogateBoundaryModeler::ClearGapModelPart2D()
{
    const std::string output_name =
        mGapParameters["output_model_part_name"].GetString();
    KRATOS_ERROR_IF_NOT(mpGapModel->HasModelPart(output_name))
        << "GapSurrogateBoundaryModeler::SetupModelPart must be called before "
        << "PrepareGeometryModel." << std::endl;
    ModelPart& r_output = mpGapModel->GetModelPart(output_name);

    const std::string gap_elements_name =
        mGapParameters["gap_element_sub_model_part_name"].GetString();
    if (r_output.HasSubModelPart(gap_elements_name)) {
        ModelPart& r_gap_elements = r_output.GetSubModelPart(gap_elements_name);
        std::vector<std::size_t> element_ids;
        element_ids.reserve(r_gap_elements.NumberOfElements());
        for (const auto& r_element : r_gap_elements.Elements()) {
            element_ids.push_back(r_element.Id());
        }
        for (const std::size_t element_id : element_ids) {
            r_output.RemoveElementFromAllLevels(element_id);
        }
        r_gap_elements.Nodes().clear();
    }

    const std::string gap_interfaces_name =
        mGapParameters["gap_interface_sub_model_part_name"].GetString();
    if (r_output.HasSubModelPart(gap_interfaces_name)) {
        ModelPart& r_gap_interfaces =
            r_output.GetSubModelPart(gap_interfaces_name);
        std::vector<std::size_t> condition_ids;
        condition_ids.reserve(r_gap_interfaces.NumberOfConditions());
        for (const auto& r_condition : r_gap_interfaces.Conditions()) {
            condition_ids.push_back(r_condition.Id());
        }
        for (const std::size_t condition_id : condition_ids) {
            r_output.RemoveConditionFromAllLevels(condition_id);
        }
        r_gap_interfaces.Nodes().clear();
    }
}

void GapSurrogateBoundaryModeler::SetupGapModelPart2D()
{
    ModelPart& r_output = mpGapModel->GetModelPart(mGapParameters["output_model_part_name"].GetString());
    ModelPart& r_surrogate = r_output.GetSubModelPart(mGapParameters["surrogate_sub_model_part_name"].GetString());
    const auto gap_name = mGapParameters["gap_element_sub_model_part_name"].GetString();
    const auto interface_name = mGapParameters["gap_interface_sub_model_part_name"].GetString();
    ModelPart& r_gap_elements = r_output.HasSubModelPart(gap_name) ? r_output.GetSubModelPart(gap_name) : r_output.CreateSubModelPart(gap_name);
    ModelPart& r_gap_interfaces = r_output.HasSubModelPart(interface_name) ? r_output.GetSubModelPart(interface_name) : r_output.CreateSubModelPart(interface_name);
    const auto lower = mGapParameters["lower_point"].GetVector();
    const auto upper = mGapParameters["upper_point"].GetVector();
    const double tolerance = 1.0e-12 * std::max({upper[0]-lower[0], upper[1]-lower[1], 1.0});
    GapConstructionUtility::Build(r_output, r_surrogate, r_gap_elements, r_gap_interfaces,
        mGapParameters["gap_element_name"].GetString(),
        mGapParameters["gap_interface_condition_name"].GetString(), tolerance,
        mGapParameters["echo_level"].GetInt());
}

} // namespace Kratos
